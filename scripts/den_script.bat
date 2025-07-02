@echo off

MOVE INITIAL.txt simulation\INITIAL.txt
REM Run the Python script to generate the adjacency matrix
python matrixGeneration\adjacencymat.py
echo Adjacency matrix generated.
REM Move the adjacency_matrix.txt file to the simulation folder
MOVE adjacency_matrix.txt simulation\adjacency_matrix.txt

REM Check if the Python script ran successfully
if %errorlevel% neq 0 (
    echo Python script failed.
    exit /b %errorlevel%
)

REM Compile the Fortran program (if not already compiled)
if not exist simulation\simulation.exe (
    cd simulation\
    gfortran -o simulation.exe sync.f90 r1279.f90 ran2.f
    cd ..
)


REM Run the Fortran program
cd simulation\
simulation.exe
cd ..

REM Check if the Fortran program ran successfully
if %errorlevel% neq 0 (
    echo Fortran program failed.
    exit /b %errorlevel%
)
REM Move output.txt and adjacency_matrix.txt to correlationCalc folder
REM MOVE simulation\output.txt correlationCalc\output.txt

