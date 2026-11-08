@echo off
setlocal
cd /d "%~dp0"
REM Validate every candidate; the first python on PATH may belong to an
REM unrelated installation. An explicit override must work or report why.
set "INSILICO_PY="
set "INSILICO_EXIT=0"
if defined INSILICO_PYTHON (
    call :try_python "%INSILICO_PYTHON%"
    if not defined INSILICO_PY goto custom_error
    goto ready
)
if defined CONDA_PREFIX call :try_python "%CONDA_PREFIX%\python.exe"
if defined INSILICO_PY goto ready
for %%P in (
    "%~dp0python.exe"
    "%USERPROFILE%\anaconda3\envs\insilico\python.exe"
    "%USERPROFILE%\miniconda3\envs\insilico\python.exe"
    "%ProgramData%\anaconda3\envs\insilico\python.exe"
    "%ProgramData%\miniconda3\envs\insilico\python.exe"
    "%USERPROFILE%\anaconda3\python.exe"
    "%USERPROFILE%\miniconda3\python.exe"
    "%ProgramData%\anaconda3\python.exe"
    "%ProgramData%\miniconda3\python.exe"
) do call :try_python %%P
if defined INSILICO_PY goto ready
echo A suitable Python environment was not found.
echo From Anaconda Prompt in this folder, run:
echo     conda env create -f environment.yml
echo     conda activate insilico
echo     python app.py
echo Or set INSILICO_PYTHON to the full path of a compatible python.exe.
echo See INSTALL.md for AutoDock Vina on Windows.
pause
exit /b 1
:custom_error
echo The Python selected by INSILICO_PYTHON could not run this app:
echo     %INSILICO_PYTHON%
echo Check that the file exists and has Python 3.11 or newer, NiceGUI
echo 3.18 or newer in version 3, RDKit, OpenMM, PDBFixer, Meeko, NumPy and SciPy.
echo Clear INSILICO_PYTHON to let the launcher search automatically.
pause
exit /b 1
:ready
echo In-Silico Docking Studio
echo Using: %INSILICO_PY%
if "%INSILICO_PORT%"=="" set "INSILICO_PORT=8111"
set "INSILICO_URL=http://127.0.0.1:%INSILICO_PORT%"
if /I "%~1"=="workshop" set "INSILICO_URL=%INSILICO_URL%/workshop"
if /I "%~1"=="check" goto check
echo Opening %INSILICO_URL%
echo Leave this window open while you work.
start "" "%INSILICO_URL%"
"%INSILICO_PY%" app.py
set "INSILICO_EXIT=%ERRORLEVEL%"
goto end
:check
"%INSILICO_PY%" -m core.selftest --offline
set "INSILICO_EXIT=%ERRORLEVEL%"
pause
:end
endlocal & exit /b %INSILICO_EXIT%

:try_python
if defined INSILICO_PY exit /b 0
if not exist "%~1" exit /b 1
"%~1" -c "import sys,nicegui,rdkit,openmm,pdbfixer,meeko,numpy,scipy; from importlib.metadata import version; v=tuple(int(n) for n in version('nicegui').split('.')[:2]); assert sys.version_info >= (3,11) and (3,18) <= v < (4,0)" >nul 2>nul
if errorlevel 1 exit /b 1
set "INSILICO_PY=%~1"
exit /b 0
