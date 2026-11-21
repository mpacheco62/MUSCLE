@echo off
rem SPDX-License-Identifier: GPL-3.0-or-later
rem Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
rem
rem Opens VS Code on this repository with the build environment of MUSCLE on Windows (README, Getting started):
rem Intel oneAPI (setvars.bat, when installed) and then the conda environment of environment.yml, so that CMake Tools
rem finds ifx, icx, gfortran, CMake, Ninja and fypp with any preset. The conda environment goes last: its CMake comes
rem before the one of Visual Studio that setvars.bat puts on the PATH.
rem   tools\vscode.cmd                 (or double click)   conda environment "muscle"
rem   tools\vscode.cmd myenv                               another conda environment
rem   tools\vscode.cmd --check [env]                       only show which tools it finds, without opening VS Code
setlocal
set "MODE=open"
if /i "%~1"=="--check" (set "MODE=check" & shift)
set "ENV_NAME=%~1"
if "%ENV_NAME%"=="" set "ENV_NAME=muscle"
for %%I in ("%~dp0..") do set "REPO=%%~fI"

rem 1. Intel oneAPI (optional: the gfortran presets do not need it)
set "ONEAPI=%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat"
if exist "%ONEAPI%" (
    echo Loading Intel oneAPI ...
    call "%ONEAPI%" intel64 > nul 2>&1
    call :msvc
) else (
    echo Intel oneAPI not found: only the gfortran presets will work.
)

rem 2. The conda environment (the conda of the console, or a usual installation folder)
set "CONDA_BAT="
if defined CONDA_EXE for %%D in ("%CONDA_EXE%\..\..") do if exist "%%~fD\condabin\conda.bat" set "CONDA_BAT=%%~fD\condabin\conda.bat"
if not defined CONDA_BAT for %%R in ("%USERPROFILE%" "%LOCALAPPDATA%" "%ProgramData%") do (
    for %%N in (miniforge3 miniconda3 anaconda3) do (
        if not defined CONDA_BAT if exist "%%~R\%%N\condabin\conda.bat" set "CONDA_BAT=%%~R\%%N\condabin\conda.bat"
    )
)
if not defined CONDA_BAT (
    echo ERROR: conda not found. Install Miniforge ^(README, Getting started, step 1^).
    goto :fail
)
call "%CONDA_BAT%" activate %ENV_NAME%
if errorlevel 1 (
    echo ERROR: could not activate the conda environment "%ENV_NAME%". Create it with:
    echo     conda env create -f environment.yml
    goto :fail
)
echo Conda environment: %CONDA_PREFIX%

if "%MODE%"=="check" (
    for %%T in (gfortran ifx icx link cmake ninja fypp ford code) do call :show %%T
    exit /b 0
)

rem 3. VS Code: a window that is already open keeps its old environment
echo Opening VS Code on %REPO%
echo (if VS Code was already open, close every window first and run this again)
call code "%REPO%"
exit /b 0

:fail
pause
exit /b 1

rem ifx and icx link with the Microsoft linker. setvars.bat loads it (vcvars64.bat) only when it recognizes the
rem Visual Studio installation; otherwise ifx fails with "could not find 'link.exe'" and icx with "cannot open file
rem 'libcmt.lib'". Then load it from the Visual Studio (or Build Tools) with the C++ tools that vswhere finds.
:msvc
where link > nul 2>&1 && exit /b 0
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
set "VSDIR="
if exist "%VSWHERE%" for /f "usebackq delims=" %%V in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSDIR=%%V"
rem No ( ) block here: the path usually has "(x86)", whose ")" would close the block
if not defined VSDIR goto :msvc_check
if not exist "%VSDIR%\VC\Auxiliary\Build\vcvars64.bat" goto :msvc_check
echo Loading the Microsoft C++ tools of "%VSDIR%" ...
call "%VSDIR%\VC\Auxiliary\Build\vcvars64.bat" > nul 2>&1
:msvc_check
where link > nul 2>&1 && exit /b 0
echo WARNING: no Microsoft linker ^(link.exe^): ifx and icx cannot link. Install the C++ tools of Visual Studio
echo          ^(README, Getting started, Intel ifx: Build Tools with the workload "Desktop development with C++"^).
exit /b 0

rem The first match of a tool on the PATH
:show
for /f "delims=" %%P in ('where %1 2^>nul') do (
    echo   %1: %%P
    exit /b 0
)
echo   %1: not found
exit /b 0
