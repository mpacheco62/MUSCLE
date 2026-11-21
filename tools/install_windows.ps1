# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Installs on Windows what MUSCLE needs to build (README, Getting started): Git, Miniforge and the conda environment of
# environment.yml (GNU Fortran, CMake, Ninja, OpenBLAS, fypp, FORD); optionally VS Code with its extensions and, for the
# ifx-* presets, the C++ tools of Visual Studio (Build Tools), the Intel Fortran Compiler and oneMKL. What is already
# installed is skipped. The Intel C/C++ compiler (icx) is not in winget: the script says how to install it by hand.
#
#   powershell -ExecutionPolicy Bypass -File tools\install_windows.ps1 [-VSCode] [-Ifx] [-EnvName muscle] [-DryRun]
#
# -DryRun only shows what it would install. The installers of Visual Studio and Intel ask for administrator rights.
param(
    [switch]$VSCode,          # Visual Studio Code and the extensions CMake Tools and Modern Fortran
    [switch]$Ifx,             # Build Tools (C++), Intel Fortran Compiler and oneMKL; then icx by hand
    [string]$EnvName = 'muscle',
    [switch]$DryRun
)
$ErrorActionPreference = 'Stop'
$Repo = Split-Path -Parent $PSScriptRoot
$OneApi = Join-Path ${env:ProgramFiles(x86)} 'Intel\oneAPI'
$VsWhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
$Pending = New-Object System.Collections.Generic.List[string]

function Say([string]$Text) { Write-Host "== $Text" -ForegroundColor Cyan }

# The PATH of the registry, so that what winget just installed is found in this same console
function Update-SessionPath {
    $machine = [Environment]::GetEnvironmentVariable('Path', 'Machine')
    $user = [Environment]::GetEnvironmentVariable('Path', 'User')
    $env:Path = "$machine;$user"
}

# Runs a program without showing its output and returns its exit code. Windows PowerShell 5.1 turns every stderr line
# of a program whose output is redirected into an error, which stops the script with $ErrorActionPreference = 'Stop'
# (code --install-extension prints a Node warning there): only the exit code decides here.
function Invoke-Quiet([string]$Exe, [string[]]$ArgList) {
    $previous = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try { & $Exe @ArgList 2>&1 | Out-Null } finally { $ErrorActionPreference = $previous }
    return $LASTEXITCODE
}

function Test-WingetPackage([string]$Id) {
    return ((Invoke-Quiet 'winget' @('list', '--id', $Id, '-e', '--accept-source-agreements')) -eq 0)
}

function Install-WingetPackage([string]$Id, [string]$Name, [string[]]$Extra = @()) {
    if (Test-WingetPackage $Id) { Write-Host "   $Name ($Id): already installed"; return }
    $cmd = @('install', '--id', $Id, '-e', '--accept-package-agreements', '--accept-source-agreements') + $Extra
    if ($DryRun) { Write-Host "   would run: winget $($cmd -join ' ')"; return }
    Write-Host "   installing $Name ($Id) ..."
    & winget @cmd
    if ($LASTEXITCODE -ne 0) { throw "winget could not install $Id (exit code $LASTEXITCODE)" }
    Update-SessionPath
}

# conda.bat of the first usual installation (Miniforge first)
function Find-Conda {
    foreach ($root in @($env:USERPROFILE, $env:LOCALAPPDATA, $env:ProgramData)) {
        foreach ($name in @('miniforge3', 'miniconda3', 'anaconda3')) {
            $bat = Join-Path $root "$name\condabin\conda.bat"
            if (Test-Path $bat) { return $bat }
        }
    }
    return $null
}

# 0. winget
Say 'winget'
if (-not (Get-Command winget -ErrorAction SilentlyContinue)) {
    Write-Host '   winget is missing. Install "App Installer" from the Microsoft Store or, in a PowerShell as administrator:'
    Write-Host '     Install-PackageProvider -Name NuGet -Force'
    Write-Host '     Install-Module -Name Microsoft.WinGet.Client -Force -Repository PSGallery'
    Write-Host '     Repair-WinGetPackageManager -AllUsers'
    Write-Host '   then open a new console and run this script again.'
    exit 1
}
Write-Host "   $(winget --version)"

# 1. Git and Miniforge
Say 'Git and Miniforge'
if (Get-Command git -ErrorAction SilentlyContinue) { Write-Host '   Git: already installed' }
else { Install-WingetPackage 'Git.Git' 'Git' }
$conda = Find-Conda
if ($conda) { Write-Host "   conda: already installed ($conda)" }
else { Install-WingetPackage 'CondaForge.Miniforge3' 'Miniforge'; $conda = Find-Conda }

# 2. The conda environment of environment.yml
Say "Conda environment '$EnvName'"
if (-not $conda) {
    if ($DryRun) { Write-Host "   would run: conda env create -f environment.yml -n $EnvName" }
    else { throw 'conda not found after installing Miniforge: open a new console and run this script again' }
} else {
    $envDir = Join-Path (Split-Path -Parent (Split-Path -Parent $conda)) "envs\$EnvName"
    if (Test-Path $envDir) {
        Write-Host "   already exists ($envDir); to update it: conda env update -n $EnvName -f environment.yml --prune"
    } elseif ($DryRun) {
        Write-Host "   would run: conda env create -f environment.yml -n $EnvName"
    } else {
        & $conda env create -f (Join-Path $Repo 'environment.yml') -n $EnvName
        if ($LASTEXITCODE -ne 0) {
            throw ('conda env create failed. With Miniconda, accept the terms of the "defaults" channel first ' +
                   '(conda tos accept) or remove it (conda config --remove channels defaults)')
        }
    }
    # conda activate in PowerShell (adds a block to the profile; does nothing when it is already there)
    if ($DryRun) { Write-Host '   would run: conda init powershell' }
    else { Invoke-Quiet $conda @('init', 'powershell') | Out-Null }
    $policy = Get-ExecutionPolicy
    if ($policy -eq 'Restricted' -or $policy -eq 'AllSigned') {
        $Pending.Add('PowerShell may refuse to run the conda profile: Set-ExecutionPolicy -Scope CurrentUser RemoteSigned')
    }
}

# 3. VS Code (optional)
if ($VSCode) {
    Say 'Visual Studio Code'
    Install-WingetPackage 'Microsoft.VisualStudioCode' 'VS Code'
    $code = Join-Path $env:LOCALAPPDATA 'Programs\Microsoft VS Code\bin\code.cmd'
    if (-not (Test-Path $code)) { $code = (Get-Command code -ErrorAction SilentlyContinue).Source }
    foreach ($ext in @('ms-vscode.cmake-tools', 'fortran-lang.linter-gfortran')) {
        if ($DryRun -or -not $code) { Write-Host "   would run: code --install-extension $ext"; continue }
        if ((Invoke-Quiet $code @('--install-extension', $ext, '--force')) -eq 0) { Write-Host "   extension $ext" }
        else { $Pending.Add("Install the VS Code extension $ext (code --install-extension $ext)") }
    }
}

# 4. ifx (optional): C++ tools of Visual Studio, Intel Fortran Compiler, oneMKL; icx by hand
if ($Ifx) {
    Say 'C++ tools of Visual Studio (link.exe and the Windows SDK, needed by ifx)'
    $vcTools = 'Microsoft.VisualStudio.Component.VC.Tools.x86.x64'
    $withVc = $null
    if (Test-Path $VsWhere) { $withVc = & $VsWhere -latest -products * -requires $vcTools -property installationPath }
    if ($withVc) {
        Write-Host "   already installed ($withVc)"
    } else {
        $buildTools = $null
        if (Test-Path $VsWhere) {
            $buildTools = & $VsWhere -latest -products Microsoft.VisualStudio.Product.BuildTools -property installationPath
        }
        $workload = @('--add', 'Microsoft.VisualStudio.Workload.VCTools', '--includeRecommended', '--passive', '--wait')
        if ($buildTools) {
            # Build Tools without the C++ workload: add it to that installation
            $setup = Join-Path (Split-Path -Parent $VsWhere) 'setup.exe'
            $cmd = @('modify', '--installPath', $buildTools) + $workload
            if ($DryRun) { Write-Host "   would run: `"$setup`" $($cmd -join ' ')" }
            else {
                Write-Host "   adding the C++ workload to $buildTools ..."
                $p = Start-Process -FilePath $setup -ArgumentList $cmd -Wait -PassThru
                if ($p.ExitCode -ne 0 -and $p.ExitCode -ne 3010) { throw "Visual Studio Installer failed ($($p.ExitCode))" }
            }
        } else {
            Install-WingetPackage 'Microsoft.VisualStudio.2022.BuildTools' 'Build Tools 2022' @('--override', ($workload -join ' '))
        }
    }

    Say 'Intel Fortran Compiler and oneMKL'
    if (Test-Path (Join-Path $OneApi 'compiler\latest\bin\ifx.exe')) { Write-Host '   ifx: already installed' }
    else { Install-WingetPackage 'Intel.FortranCompiler' 'Intel Fortran Compiler' }
    if (Test-Path (Join-Path $OneApi 'mkl\latest')) { Write-Host '   oneMKL: already installed' }
    else { Install-WingetPackage 'Intel.oneMKL' 'Intel oneMKL' }

    Say 'Intel C/C++ compiler (icx)'
    if (Test-Path (Join-Path $OneApi 'compiler\latest\bin\icx.exe')) { Write-Host '   icx: already installed' }
    else {
        Write-Host '   not installed, and it is not in winget. Install it by hand (same version as ifx):'
        Write-Host '     https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler-download.html'
        $Pending.Add('Install the Intel oneAPI DPC++/C++ Compiler (icx), needed by the ifx-* presets (link above)')
    }
}

# 5. What the build environment finds (tools\vscode.cmd --check), and what is left
Say 'Check'
if ($DryRun) { Write-Host '   would run: tools\vscode.cmd --check' }
else { & cmd /c "`"$(Join-Path $Repo 'tools\vscode.cmd')`" --check $EnvName" }
if ($Pending.Count -gt 0) {
    Say 'Pending'
    foreach ($item in $Pending) { Write-Host "   - $item" -ForegroundColor Yellow }
}
Write-Host ''
Write-Host 'Done. Open a new console (so that it sees the new PATH) and continue with the README, Getting started, step 3;'
Write-Host 'or open VS Code with tools\vscode.cmd.'
