# MUSCLE: Mechanics User-defined Solids Constitutive Library & Extensions

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: Fortran](https://img.shields.io/badge/Language-Fortran_2008-orange.svg)](https://gcc.gnu.org/fortran/)

**MUSCLE** is a high-performance Fortran 2008 library for the development, testing and implementation of
constitutive material models in Computational Solid Mechanics (elasticity, plasticity, viscoplasticity,
viscoelasticity, finite strain kinematics).

The sources are [FYPP](https://fypp.readthedocs.io) templates: every tensor type, model and material is
generated in single (`sp`) and double (`dp`) precision from one source, with static dispatch only (no
polymorphism, `pure` procedures), so that the same code runs on CPUs (OpenMP) and can be compiled for GPUs
(OpenMP offload). Constitutive models are combined into **declared materials**: you list the combinations
you need in a small Python registry and the build generates one Fortran module per material, with a point
API for the user subroutines of FE codes (Abaqus UMAT, ANSYS USERMAT, LS-DYNA umat41) and a batch API over
arrays of points.

---

## Key Features

### 1. Tensor Engine (`muscle_tensors`)
Voigt-notation tensor algebra with manually unrolled loops, in 3D and 2D, `sp` and `dp`:
*   **2nd order tensors**, symmetric and general (`ten_3D2Osym_dp`, `ten_2D2O_sp`, ...).
*   **4th order tensors** with minor, major and full symmetries (`ten_3D4O3sym_dp`, ...).
*   **Identities**, plain and scaled (`iden_2O`, `iden_4O4T`, `iden_2OS_dp`, ...).
*   **Overloaded operators:** `+ - * / =`, double contraction (`.ddot.`), dyadic products (`.tdot.`,
    `.tdotsym.`), `.transform.`, deviatoric part (`.dev.`) and inverse (`.inv.`).
*   Names without kind suffix (`ten_3D2Osym`, ...) are the `dp` types.

### 2. Model families (free `pure` functions over concrete parameter types)
*   **Elasticity:** linear isotropic; frozen viscoelastic elasticity driven by a viscoelastic model.
*   **Viscoelasticity:** generalized Maxwell and Kelvin-Voigt.
*   **Yield criteria:** von Mises, Hill48, Drucker-Prager, CPB06, CPB06ex2 (Cazacu-Plunkett-Barlat) and
    Yld2004-18p, with analytical first and second derivatives.
*   **Isotropic hardening:** bilinear, Ludwik, Swift, modified Voce.
*   **Overstress viscoplasticity:** Perzyna, Chaboche.
*   **Viscoplastic flow stress laws:** Johnson-Cook (and full Johnson-Cook), Rusinek-Klepaczko (standard and
    modified), Nemat-Nasser-Li, Voyiadjis-Almasri.
*   **Kinematics:** small strain, total and updated Lagrangian, material and spatial logarithmic,
    corotational (strain rotated by the FE code).

### 3. Return-mapping solvers (bodies included by every material)
*   **Closest point projection** (Newton-Raphson, consistent tangent), with overstress (`closest_point_perzyna`)
    and consistency (`closest_point_time`) viscoplasticity.
*   **Monolithic closest point** (`closest_point_2`, 7x7 Newton with an own LU solver).
*   **Cutting plane** (explicit, continuum elastoplastic tangent).

Which one to declare a material with:

| Solver | For | Iterations | Tangent | When |
|--------|-----|------------|---------|------|
| `closest_point` | rate-independent plasticity | scalar Newton on the plastic multiplier \(\Delta\gamma \ge 0\), damped (x0.75) when it does not converge | consistent | **the default**: quadratic convergence |
| `closest_point_2` | rate-independent plasticity | monolithic 7x7 Newton on the plastic strain and \(\Delta\gamma\), own LU, no \(\Delta\gamma \ge 0\) | the same as `closest_point` | compatibility with the historical solver; known issue: it converges linearly with shear components (its Jacobian lacks the factor 2 of the shears), prefer `closest_point` |
| `cutting_plane` | rate-independent plasticity | explicit corrections (at most 100), no linear system | continuum, not consistent | cheap iterations; fits an explicit FE code (no tangent needed); with an implicit one the global convergence is slower |
| `closest_point_perzyna` | overstress viscoplasticity (Perzyna, Chaboche) | Newton with a line search | consistent | a material with an overstress (`os`) law |
| `closest_point_time` | consistency viscoplasticity (Johnson-Cook, ...) | Newton | consistent | a material with a flow stress (`vp`) law; known issue: Johnson-Cook with a very small overstress may not converge |

### 4. Declared materials (`muscle_mat_<name>`)
For each material of the registry, in `sp` and `dp`:
*   `mat_<name>_<k>`: the material type (one parameter component per model, `tol`, `iter_nw`).
*   `mat_init_<name>(mat, props, status)` and `MAT_<NAME>_NPROPS`: the material from a flat parameter array.
*   `muscle_point_<name>(...)`: point API over flat arrays (`strain(6)` or `F_n(9)`, `F_np1(9)`), for FE
    user subroutines.
*   `muscle_update_<name>(...)`: batch API over structure-of-arrays data (`eps(npts, 6)`), OpenMP on the
    CPU, optional index list for multiphase models.
*   `MAT_<NAME>_NSTATE`, `mat_pack_<name>`/`mat_unpack_<name>`, `mat_rotate_hsv_<name>`: the history layout.

---

## Requirements

*   **Compiler:** a Fortran 2008 compiler with OpenMP (the GPU build needs OpenMP 5.1 offload). Older
    versions than the tested ones below may work but are not verified.
*   **Build:** CMake 3.20+ (Ninja recommended) and **FYPP 3.1+** (`pip install fypp`).
*   **Python 3** (FYPP and the material registry).
*   **BLAS & LAPACK** (or Intel MKL), required by CMake.
*   Optional: [FORD](https://forddocs.readthedocs.io) for the documentation.

The conda environment of `environment.yml` (conda-forge) has all of them, on Linux and on Windows: GNU Fortran,
CMake, Ninja, BLAS/LAPACK, Python, FYPP and FORD (`conda env create -f environment.yml`, then
`conda activate muscle`). Intel `ifx` comes from oneAPI, outside the environment. The environment takes its
packages from conda-forge only; on a Miniconda whose configuration still lists the Anaconda `defaults` channel,
conda asks first to accept that channel's terms of service (`conda tos accept`) or to remove it
(`conda config --remove channels defaults`); Miniforge needs neither.

Tested toolchain (Linux; on Windows, every test passes with GNU Fortran 16.2 of conda-forge and with `ifx` 2026.2):

| Tool | Version | Scope |
|------|---------|-------|
| GNU Fortran | 15.2 | CPU build and every test (`sp;dp`, Debug); NVPTX and AMDGCN offload builds compile and link (not run on NVIDIA/AMD hardware) |
| Intel `ifx` (oneAPI) | 2026.2 | SPIR64 offload build; target tests run on an Intel UHD 770 GPU (some cases excluded by IGC bugs, `tools/gpu_skip_tests.sh`) |
| CMake | 4.2 (3.20+ required) | |
| Ninja | 1.13 | |
| FYPP | 3.1 and 3.3 | |
| Python | 3.10 and 3.14 | |

Formatted output of the tensors (UDTIO) is disabled by default with GNU Fortran < 16.3, whose UDTIO is
not thread-safe under OpenMP (`MUSCLE_UDTIO=AUTO`).

---

## Getting started: from a fresh clone to your own program

The commands work the same in a Linux shell and in Windows PowerShell (keep the quotes around the `-D` options that
contain `$PWD` or `$HOME`: PowerShell does not expand a variable inside an unquoted `-D...` argument).

**1. Install the tools.** Git and [Miniforge](https://conda-forge.org/download/) (conda with the conda-forge channel;
Miniconda also works, see Requirements). The compiler, CMake and the rest come in step 2.

On Windows there are two ways: the script `tools\install_windows.ps1` (installs everything in one command), or by hand
(the commands below, then steps 2 and 3, and *Optional: Intel ifx*). Both end in the same place.

> **Windows with the script.** Needs `winget` (see below if it is missing), Git (`winget install --id Git.Git -e`) and
> the repository cloned (step 2, the `git clone`). From the folder of the repository:
>
> ```powershell
> powershell -ExecutionPolicy Bypass -File tools\install_windows.ps1 -VSCode -Ifx
> ```
>
> It installs with `winget` what is missing and skips what is already there: Miniforge and the conda environment
> `muscle` (step 2), and optionally VS Code with its extensions (`-VSCode`) and what the `ifx-*` presets need (`-Ifx`:
> Build Tools with the C++ tools, Intel Fortran Compiler, oneMKL). The Intel C/C++ compiler (`icx`) is not in
> `winget`: the script lists it as pending, with the page to download it. At the end it runs `tools\vscode.cmd --check`
> and lists what is left to do. `-DryRun` only shows what it would do; `-EnvName <name>` uses another environment. The
> installers of Visual Studio and Intel ask for administrator rights. Then open a new console and continue with step 3
> (or open VS Code with `tools\vscode.cmd`, *Working in VS Code*).

*Windows by hand* (PowerShell). The commands use `winget`, which comes with Windows 10 (1809 or later) and Windows 11 as part of
*App Installer*. Check it with `winget --version`; if it is missing, install *App Installer* from the Microsoft Store, or
(Windows Server, or without the Store) from a PowerShell opened as administrator:

```powershell
Install-PackageProvider -Name NuGet -Force
Install-Module -Name Microsoft.WinGet.Client -Force -Repository PSGallery
Repair-WinGetPackageManager -AllUsers          # installs winget and its dependencies; then open a new console
```

Without `winget`, the installers of [Git](https://git-scm.com/downloads/win) and
[Miniforge](https://conda-forge.org/download/) work the same. Then:

```powershell
winget install --id Git.Git -e
winget install --id CondaForge.Miniforge3 -e
```

Then open the *Miniforge Prompt* (start menu) and run `conda init powershell` once, so that `conda activate` works in
PowerShell; close and open PowerShell again (if it refuses to run the profile script:
`Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`).

*Linux* (Debian/Ubuntu; other distributions: their package manager for Git):

```bash
sudo apt install git
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh      # accept the license, and "yes" to initialize conda; reopen the shell
```

*Optional:*
- Visual Studio Code (Windows: `winget install --id Microsoft.VisualStudioCode -e`; Linux: the package of
  [code.visualstudio.com](https://code.visualstudio.com/download)) with the extensions *CMake Tools* and *Modern Fortran*:
  `code --install-extension ms-vscode.cmake-tools` and `code --install-extension fortran-lang.linter-gfortran`.
- Intel `ifx`: see below.

Check: `git --version` and `conda --version` answer in a new console.

*Optional: Intel `ifx`* (the `ifx-*` presets). GNU Fortran from step 2 is enough to build and use MUSCLE.

- **Windows.** `ifx` compiles, but links with the Microsoft linker (`link.exe`) and the Windows SDK, which come with
  Visual Studio; without them the oneAPI installer warns that it finds no Visual Studio and `ifx` cannot link.
  1. The C++ tools of Visual Studio (tested with Visual Studio 2022; the free *Build Tools* are enough, no IDE):

     ```powershell
     winget install --id Microsoft.VisualStudio.2022.BuildTools -e `
         --override "--passive --wait --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
     ```

     (or, in the *Visual Studio Installer* of an existing Visual Studio, the workload *Desktop development with C++*).
  2. **Two** Intel compilers, each from its download page (Windows, *Online* or *Offline* installer): the
     [Intel Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler-download.html)
     (`ifx`) and the
     [Intel oneAPI DPC++/C++ Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler-download.html)
     (`icx`, the C compiler of the `ifx-*` presets; the Fortran installer does not include it, and without it CMake
     stops with `CMAKE_C_COMPILER: icx is not a full path and was not found in the PATH`). Both install in
     `C:\Program Files (x86)\Intel\oneAPI`; check that `compiler\latest\bin` has `ifx.exe` and `icx.exe`. If you
     installed them before Visual Studio, there is no need to reinstall them: `setvars.bat` finds Visual Studio when it
     runs.
  3. Intel oneMKL, the BLAS/LAPACK of `ifx` on Windows (`winget install --id Intel.oneMKL -e`, or its download page).
     The OpenBLAS of the conda environment does not work with `ifx` on Windows: it exports the names of GNU Fortran
     (`sgemm_`) and `ifx` calls `SGEMM`, so CMake stops with `Could NOT find BLAS (missing: BLAS_LIBRARIES)`. On Linux
     both compilers use the same names and the OpenBLAS of the environment is enough.
  4. Every time you use `ifx`, oneAPI must be loaded, and before the conda environment (`setvars.bat` puts the CMake of
     Visual Studio on the `PATH`; the one of the environment must come first). In VS Code, `tools\vscode.cmd` does it
     (see *Working in VS Code* below); in a terminal, load oneAPI into a new PowerShell and activate the environment
     there:

     ```powershell
     cmd /c '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'
     conda activate muscle
     ifx /QV                          # Intel(R) Fortran Compiler ... Version 2026.x
     icx /QV                          # Intel(R) oneAPI DPC++/C++ Compiler ... Version 2026.x
     cmake --preset ifx-release       # and the rest of step 3 with ifx-release
     ```

     If CMake then fails with `ifx: error #10037: could not find 'link.exe'` or `icx: ... cannot open file
     'libcmt.lib'`, `setvars.bat` did not recognize the Visual Studio installation (it happens with the Build Tools
     alone): load its C++ tools by hand before `conda activate`, with
     `cmd /c '"C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" && powershell'`
     (`...\2022\Community\...` for Visual Studio Community), and check with `where.exe link`. `tools\vscode.cmd` does
     it by itself.
- **Linux.** The same two compilers (no Visual Studio needed), from the download pages or from Intel's APT repository
  (packages `intel-oneapi-compiler-fortran` and `intel-oneapi-compiler-dpcpp-cpp`); every time:
  `source /opt/intel/oneapi/setvars.sh`, then the `ifx-*` presets.

**2. Get the code and create the environment** (GNU Fortran, CMake, Ninja, OpenBLAS, FYPP, FORD; once per machine):

```bash
git clone https://github.com/mpacheco62/muscle.git
cd muscle
conda env create -f environment.yml
conda activate muscle               # in every new console
```

**3. Build the library and run the tests:**

```bash
cmake --preset gfortran-release     # configure (build/gfortran-release)
cmake --build --preset gfortran-release
ctest --preset gfortran-release     # every test should pass
```

*Working in VS Code.* The extension CMake Tools uses the tools of the environment VS Code was started from, and
VS Code started from the start menu has neither the conda environment nor oneAPI (errors such as `icx is not a full
path and was not found in the PATH`, `fypp not found`, or the CMake of Visual Studio instead of the one of the
environment). So:

- **Windows:** close every VS Code window and open it with `tools\vscode.cmd` (double click in the Explorer, or run it
  from a console). It loads Intel oneAPI when it is installed, then the conda environment `muscle`, and opens VS Code
  on the repository. `tools\vscode.cmd --check` only lists the tools it finds (`gfortran`, `ifx`, `icx`, `link`,
  `cmake`, `ninja`, `fypp`, ...), to see what is missing; `tools\vscode.cmd <env>` uses another conda environment.
- **Linux:** `conda activate muscle` (and `source /opt/intel/oneapi/setvars.sh` for `ifx`), then `code .`.

Then *CMake: Select Configure Preset* (for example *gfortran Debug* or *ifx Release*), `F7` builds, and
*CMake: Run Tests* runs the tests. After changing the way VS Code is opened, *CMake: Delete Cache and Reconfigure*.

**4. (Optional) Choose your materials.** The library contains one module per material of a registry; the default one,
`tools/fypp/muscle_materials.py`, has a test material per combination of models, `cp_vm_swift` among them. To build
yours, write a registry (see [Declaring materials](#declaring-materials)) and configure with it:

```bash
cmake --preset gfortran-release "-DMUSCLE_MATERIALS=$PWD/my_materials.py"
cmake --build --preset gfortran-release
```

**5. Install it** where your projects will find it (any folder; the library goes to `lib/`, the Fortran modules to
`include/`, the CMake package to `lib/cmake/muscle/`):

```bash
cmake --install build/gfortran-release --prefix "$HOME/muscle-install"
```

**6. Use it from your own CMake project.** Two files in a new folder, `myapp/`: the `main.f90` of the
[Quick Usage Example](#quick-usage-example) and this `CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.20)
project(myapp LANGUAGES Fortran)
find_package(muscle REQUIRED)        # the library, its modules, BLAS/LAPACK and OpenMP
add_executable(myapp main.f90)
target_link_libraries(myapp PRIVATE muscle)
```

```bash
cd myapp
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release "-DCMAKE_PREFIX_PATH=$HOME/muscle-install"
cmake --build build
./build/myapp                       # sigma_xx = 848.83...  iterations: 4
```

The example uses `cp_vm_swift`, from the default registry; with your own registry (step 4), use the names of your
materials. Compile your program with the same compiler that built MUSCLE (the `.mod` files belong to that compiler): with
the `ifx-*` presets, add `-DCMAKE_Fortran_COMPILER=ifx`.

---

## Building

```bash
git clone https://github.com/mpacheco62/muscle.git
cd muscle
conda env create -f environment.yml     # once
conda activate muscle

# Configure, build the library and the tests, run the tests (build/gfortran-debug)
cmake --preset gfortran-debug
cmake --build --preset gfortran-debug
ctest --preset gfortran-debug
```

`CMakePresets.json` has `gfortran-debug`, `gfortran-release`, `ifx-debug` and `ifx-release`, each in
`build/<preset>`. The `ifx` presets need the oneAPI environment loaded first (Linux:
`source /opt/intel/oneapi/setvars.sh`; Windows: the *Intel oneAPI command prompt*, or `setvars.bat`), with the
conda environment active too (for FYPP).

Without presets, with the tools of your choice on the `PATH`:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build                               # the library
cmake --build build --target build_all_tests      # the tests (not built by the default target)
ctest --test-dir build --output-on-failure
```

**Visual Studio Code** (extension CMake Tools): it lists the presets above. CMake Tools does not activate conda
environments nor load oneAPI: open VS Code with `tools\vscode.cmd` on Windows, or from a console with the environment
active on Linux (*Working in VS Code*, in Getting started). Alternatively, for the `gfortran` presets only, put the
environment in the `PATH` of CMake Tools in your `.vscode/settings.json`, for example on Windows:

```json
"cmake.environment": {
    "PATH": "${env:USERPROFILE}\\miniconda3\\envs\\muscle;${env:USERPROFILE}\\miniconda3\\envs\\muscle\\Library\\bin;${env:USERPROFILE}\\miniconda3\\envs\\muscle\\Scripts;${env:PATH}"
}
```

(on Linux: `"PATH": "${env:HOME}/miniconda3/envs/muscle/bin:${env:PATH}"`; adapt the path to your conda
installation).

Options:

| Option | Effect |
|--------|--------|
| `-DCMAKE_BUILD_TYPE=Debug\|Release\|Profile\|Profile_unittests` | Build type (`Profile_unittests`: coverage, target `coverage_report`). |
| `-DMUSCLE_MATERIALS=<path>` | Your material registry (default: the test registry `tools/fypp/muscle_materials.py`). |
| `-DMUSCLE_REAL_KINDS=dp` | Generate only double precision (default `sp;dp`). |
| `-DMUSCLE_REAL_KINDS=sp` | Generate only single precision, for GPUs without FP64 (no historical names such as `ten_3D2Osym`: use `ten_3D2Osym_sp`; only the GPU tests and the tests that do not need dp are built). |
| `-DMUSCLE_BUILD_BENCHMARKS=ON` | Build the benchmarks (`benchmarks/`, target `build_all_benchmarks`). |
| `-DMUSCLE_GPU=OFF\|NVPTX\|AMDGCN\|SPIR64` | OpenMP offload target of the kernels (default `OFF`, see below). |
| `-DMUSCLE_GPU_ARCH=<arch>` | Device architecture (NVPTX: e.g. `sm_80`; AMDGCN: default `gfx90a`; not used with SPIR64). |
| `-DMUSCLE_UDTIO=AUTO\|ON\|OFF` | Formatted output of the tensors (default `AUTO`, see Requirements). |

### GPU build (OpenMP offload)

| `MUSCLE_GPU` | Compiler | What you need installed |
|--------------|----------|-------------------------|
| `NVPTX` | GNU Fortran | GCC built with the nvptx offload backend (e.g. the distribution package `gcc-<version>-offload-nvptx`) |
| `AMDGCN` | GNU Fortran | GCC built with the amdgcn offload backend (`gcc-<version>-offload-amdgcn`) |
| `SPIR64` | Intel `ifx` | oneAPI with its environment loaded (`source /opt/intel/oneapi/setvars.sh`; Windows: `setvars.bat`) and the Intel GPU driver (Level Zero) |

```bash
# NVIDIA (gfortran with nvptx offload)
cmake -S . -B build-nvptx -G Ninja -DCMAKE_BUILD_TYPE=Release -DMUSCLE_GPU=NVPTX
# Intel GPU (ifx)
source /opt/intel/oneapi/setvars.sh
cmake -S . -B build-spir64 -G Ninja -DCMAKE_BUILD_TYPE=Release -DMUSCLE_GPU=SPIR64 \
      -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_C_COMPILER=icx
```

CMake checks the combination at configure time by linking a small program with a target region. Run
GPU programs with `OMP_TARGET_OFFLOAD=MANDATORY` to fail instead of falling back to the host; on Intel GPUs
without native double precision, `OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1` enables the
emulated FP64 on Linux. The Windows driver does not emulate FP64: there, on such a GPU, build with
`-DMUSCLE_REAL_KINDS=sp` (a device image with any double precision code does not load). On Windows, put
oneAPI's `compiler\bin` first in the `PATH` when a conda environment is active after `setvars.bat`: the
`libiomp5md.dll` of some conda environments (package `llvm-openmp`) does not work with offload. The tests run
on the GPU under a memory cap and a time limit (`tools/check.sh` on Linux; on Windows
`powershell -ExecutionPolicy Bypass -File tools\limit_run.ps1 -MemoryGB 4 -TimeoutSec 3600 -- ctest ...`).
The batch API on the device is described in [FEA integration](pages/integration.md).

Documentation: `cmake --build build --target doc` (output in `build/doc`, when CMake finds `ford`), or
`ford my_project.md` from the repository root.

---

## Declaring materials

A registry is a Python file that defines `MATERIALS` (see `tools/fypp/muscle_materials.py`, which also
documents every model and solver it covers):

```python
# my_materials.py
from muscle_meta import Material

MATERIALS = (
    Material('steel', 'closest_point', elas='linear', yield_='vonmises', hard='swift',
             doc='von Mises with Swift hardening'),
    Material('alu_rate', 'closest_point_time', elas='linear', yield_='hill48', vp='jc', hard='swift',
             doc='Hill48, Johnson-Cook, Swift hardening'),
    Material('steel_ul', 'closest_point', elas='linear', yield_='vonmises', hard='bilinear',
             kin='updated_lagrangian', tangent_rate='jaumann', doc='Updated Lagrangian'),
)
```

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DMUSCLE_MATERIALS=$PWD/my_materials.py
python3 tools/fypp/muscle_meta.py props steel my_materials.py   # order of the props of mat_init_steel
```

The registry is validated at configure time (known solver, one model of every family the solver needs,
valid names). Changing a combination of models requires rebuilding; the parameters do not.

---

## Quick Usage Example

With the material `cp_vm_swift` of the default registry (von Mises, Swift hardening, closest point):

```fortran
program main
    use muscle_kinds, only : dp
    use muscle_status, only : STATUS_OK
    use muscle_mat_cp_vm_swift, only : mat_cp_vm_swift_dp, mat_init_cp_vm_swift, muscle_point_cp_vm_swift, &
        muscle_update_cp_vm_swift, MAT_CP_VM_SWIFT_NSTATE
    implicit none

    integer, parameter :: NS = MAT_CP_VM_SWIFT_NSTATE, NPTS = 1000
    type(mat_cp_vm_swift_dp) :: mat
    real(dp) :: strain(6), stress(6), ddsdde(6, 6), hsv_n(NS), hsv_np1(NS)
    real(dp), allocatable :: eps(:, :), hsv_bn(:, :), hsv_bnp1(:, :), sig(:, :)
    integer, allocatable :: st(:), it(:)
    integer :: status, iters

    ! young, poisson, Swift k, e0, n, Newton tolerance, iteration limit (muscle_meta.py props cp_vm_swift)
    call mat_init_cp_vm_swift(mat, [210000.0_dp, 0.3_dp, 500.0_dp, 0.01_dp, 0.2_dp, 1.0e-10_dp, 200.0_dp], status)
    if (status /= STATUS_OK) error stop 'invalid parameters'

    ! One point: total strain at t_n+1, tensorial Voigt order (xx, yy, zz, xy, yz, xz), no factor 2
    strain = [0.01_dp, -0.003_dp, -0.003_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    hsv_n = 0.0_dp                       ! virgin state (plastic strain, equivalent plastic strain)
    call muscle_point_cp_vm_swift(mat, strain, 0.0_dp, hsv_n, hsv_np1, stress, status, iters, ddsdde)
    if (status == STATUS_OK) print *, 'sigma_xx =', stress(1), ' iterations:', iters

    ! A batch of points in structure-of-arrays layout (OpenMP parallel loop)
    allocate(eps(NPTS, 6), hsv_bn(NPTS, NS), hsv_bnp1(NPTS, NS), sig(NPTS, 6), st(NPTS), it(NPTS))
    eps = spread(strain, 1, NPTS)
    hsv_bn = 0.0_dp
    call muscle_update_cp_vm_swift(mat, eps, 0.0_dp, hsv_bn, hsv_bnp1, sig, st, it)
end program main
```

The same calls exist in single precision with `mat_cp_vm_swift_sp` and `real(sp)` arrays. The FE code
(or the caller) keeps the converged history `hsv_n`: the update writes `hsv_np1` and never changes
`hsv_n`, so a rejected increment is simply retried from `hsv_n`.

---

## FEA Integration

MUSCLE ships no user subroutine: each material provides the pieces an adapter needs (`mat_init_<name>`
from the material constants, `muscle_point_<name>` for the point update and the tangent `ddsdde` with the
DDSDDE shear convention, `mat_rotate_hsv_<name>` for codes that do not rotate their state variables). The
FORD page *FEA Integration* (`pages/integration.md`) has the component orders of Abaqus, ANSYS and LS-DYNA
and skeletons of UMAT, USERMAT and umat41. Coming from the former object-oriented API (`Elasticity_linear`,
`Closest_point`, `Plastic_material_history`, ...)? See `pages/migration.md`.

---

## Licensing & Commercial Use

This library is licensed under the **GNU General Public License v3.0 (GPL-3.0-or-later)**. You are free to use, modify, and redistribute this software under the terms of this license.

### Commercial Licensing Notice
> **Important Note:** If your organization wishes to integrate MUSCLE (or parts of its solver engine) into proprietary commercial software, or use it under licensing terms different from the GPLv3 (to avoid releasing your proprietary source code), please contact the author to discuss a commercial license agreement.

---

## Author & Contact

*   **Author:** Matias Pacheco-Alarcon
*   **Email:** [matias.pacheco.a@gmail.com](mailto:matias.pacheco.a@gmail.com)

*Feel free to reach out for collaborations, academic inquiries, or support with FEA integration.*
