# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Runs a command under a memory cap and a time limit on Windows: the replacement of
#   systemd-run --user --scope -p MemoryMax=<cap> -p MemorySwapMax=0  +  timeout <seconds>
# that tools/check.sh and testing.md use on Linux for every GPU run (an IGC bug once took 26 GB).
# The cap is a Job Object: it holds this PowerShell and every process it starts (ctest, the tests, the IGC JIT
# inside them), and a memory request beyond the cap fails in the process that makes it, instead of exhausting the
# machine. At the time limit the whole process tree is killed. The exit code is the one of the command
# (124 at the time limit, like timeout).
#
#   powershell -ExecutionPolicy Bypass -File tools\limit_run.ps1 [-MemoryGB 4] [-TimeoutSec 3600] `
#       [-SetEnv NAME=VALUE ...] -- ctest --test-dir build\spir64 -L target
# Everything after -- goes to the command untouched (the options of the command, such as -L, are not read here).
$ErrorActionPreference = 'Stop'
$MemoryGB = 4
$TimeoutSec = 3600
$SetEnv = @()
$Command = @()
for ($i = 0; $i -lt $args.Count; $i++) {
    switch ($args[$i]) {
        '--' { $Command = @($args[($i + 1)..($args.Count - 1)]); $i = $args.Count; break }
        '-MemoryGB' { $MemoryGB = [double]$args[++$i]; break }
        '-TimeoutSec' { $TimeoutSec = [int]$args[++$i]; break }
        '-SetEnv' { $SetEnv += $args[++$i]; break }
        default { throw "limit_run: unknown option '$($args[$i])' (the command goes after --)" }
    }
}
if ($Command.Count -eq 0) { throw 'limit_run: no command (usage: limit_run.ps1 [options] -- command args)' }

# Add-Type compiles the C# below with csc, which fails on a folder of LIB that does not exist (setvars.bat of oneAPI
# adds one): no LIB while it compiles, the command gets it back
$savedLib = $env:LIB
$env:LIB = $null
Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
public static class MuscleJob {
    [StructLayout(LayoutKind.Sequential)]
    struct Basic {
        public long PerProcessUserTimeLimit; public long PerJobUserTimeLimit; public uint LimitFlags;
        public UIntPtr MinimumWorkingSetSize; public UIntPtr MaximumWorkingSetSize; public uint ActiveProcessLimit;
        public UIntPtr Affinity; public uint PriorityClass; public uint SchedulingClass;
    }
    [StructLayout(LayoutKind.Sequential)]
    struct IoCounters {
        public ulong ReadOps, WriteOps, OtherOps, ReadBytes, WriteBytes, OtherBytes;
    }
    [StructLayout(LayoutKind.Sequential)]
    struct Extended {
        public Basic BasicLimits; public IoCounters Io; public UIntPtr ProcessMemoryLimit; public UIntPtr JobMemoryLimit;
        public UIntPtr PeakProcessMemoryUsed; public UIntPtr PeakJobMemoryUsed;
    }
    [DllImport("kernel32.dll", SetLastError = true)] static extern IntPtr CreateJobObject(IntPtr a, string name);
    [DllImport("kernel32.dll", SetLastError = true)]
    static extern bool SetInformationJobObject(IntPtr job, int cls, ref Extended info, int size);
    [DllImport("kernel32.dll", SetLastError = true)]
    static extern bool QueryInformationJobObject(IntPtr job, int cls, out Extended info, int size, IntPtr ret);
    [DllImport("kernel32.dll", SetLastError = true)] static extern bool AssignProcessToJobObject(IntPtr job, IntPtr proc);
    [DllImport("kernel32.dll")] static extern IntPtr GetCurrentProcess();
    const int ExtendedLimitInformation = 9;
    const uint LimitJobMemory = 0x200, DieOnUnhandledException = 0x400, KillOnJobClose = 0x2000;
    static IntPtr job;
    // The job, with this process in it: every process started from now on belongs to it
    public static void Enter(ulong bytes) {
        job = CreateJobObject(IntPtr.Zero, null);
        if (job == IntPtr.Zero) throw new Exception("CreateJobObject failed: " + Marshal.GetLastWin32Error());
        Extended info = new Extended();
        info.BasicLimits.LimitFlags = LimitJobMemory | DieOnUnhandledException | KillOnJobClose;
        info.JobMemoryLimit = new UIntPtr(bytes);
        if (!SetInformationJobObject(job, ExtendedLimitInformation, ref info, Marshal.SizeOf(typeof(Extended))))
            throw new Exception("SetInformationJobObject failed: " + Marshal.GetLastWin32Error());
        if (!AssignProcessToJobObject(job, GetCurrentProcess()))
            throw new Exception("AssignProcessToJobObject failed: " + Marshal.GetLastWin32Error());
    }
    // Peak committed memory of the job (it includes a request refused at the cap); 0 if the query fails
    public static ulong PeakBytes() {
        Extended info;
        if (!QueryInformationJobObject(job, ExtendedLimitInformation, out info, Marshal.SizeOf(typeof(Extended)),
                                       IntPtr.Zero)) return 0;
        return info.PeakJobMemoryUsed.ToUInt64();
    }
}
'@
$env:LIB = $savedLib

# Windows command line quoting of one argument (CommandLineToArgvW rules)
function ConvertTo-Argument([string]$Arg) {
    if ($Arg -ne '' -and $Arg -notmatch '[\s"]') { return $Arg }
    return '"' + ($Arg -replace '(\\*)"', '$1$1\"' -replace '(\\+)$', '$1$1') + '"'
}

[MuscleJob]::Enter([uint64]($MemoryGB * 1GB))
foreach ($pair in $SetEnv) {
    $name, $value = $pair -split '=', 2
    [Environment]::SetEnvironmentVariable($name, $value, 'Process')
}
$exe = (Get-Command $Command[0] -CommandType Application | Select-Object -First 1).Source
$psi = New-Object System.Diagnostics.ProcessStartInfo
$psi.FileName = $exe
$psi.Arguments = (@($Command | Select-Object -Skip 1) | ForEach-Object { ConvertTo-Argument $_ }) -join ' '
$psi.UseShellExecute = $false
$proc = [System.Diagnostics.Process]::Start($psi)
if ($proc.WaitForExit($TimeoutSec * 1000)) {
    $code = $proc.ExitCode
} else {
    Write-Host "limit_run: time limit of $TimeoutSec s reached, killing the process tree" -ForegroundColor Yellow
    # taskkill complains on stderr when a process of the tree already ended: not an error here
    $ErrorActionPreference = 'Continue'
    & taskkill.exe /T /F /PID $proc.Id 2>&1 | Out-Null
    $ErrorActionPreference = 'Stop'
    $code = 124
}
$peak = [MuscleJob]::PeakBytes() / 1GB
$inv = [Globalization.CultureInfo]::InvariantCulture
Write-Host ([string]::Format($inv, 'limit_run: exit code {0}, peak memory {1:N2} GB of {2} GB', $code, $peak, $MemoryGB))
if ($peak -gt $MemoryGB) {
    Write-Host 'limit_run: a memory request above the cap was refused (the peak counts it)' -ForegroundColor Yellow
}
exit $code
