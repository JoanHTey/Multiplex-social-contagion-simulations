@echo off
REM Run the Python script to generate the adjacency matrix
python adjacencymat.py

REM Check if the Python script ran successfully
if %errorlevel% neq 0 (
    echo Python script failed.
    exit /b %errorlevel%
)

REM Compile the Fortran program (if not already compiled)
if not exist simulation.exe (
    gfortran -o simulation.exe sync.f90 r1279.f90 ran2.f
)

REM Run the Fortran program
simulation.exe

REM Check if the Fortran program ran successfully
if %errorlevel% neq 0 (
    echo Fortran program failed.
    exit /b %errorlevel%
)

REM Run the Python script to process the results
echo Running correlation.py...
python correlation.py

REM Check if the Python script ran successfully
if %errorlevel% neq 0 (
    echo Python script correlation.py failed.
    pause
    exit /b %errorlevel%
)


echo Simulation completed successfully.