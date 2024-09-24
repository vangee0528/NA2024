@echo off
setlocal

rem Define variables
set CXX=g++
set TEX=xelatex
set SRC_DIR=.\src
set BUILD_DIR=release
set PROBLEMS=ProblemB ProblemC ProblemD ProblemE ProblemF
set TEX_SRC=report.tex
set TEX_OUT=report.pdf

@REM Create the release directory if it doesn't exist
if not exist %BUILD_DIR% (
    mkdir %BUILD_DIR%
)

@REM Main logic to execute the corresponding function based on the argument
if "%~1"=="" (
    goto :run
    goto :report
    goto :clean
    exit /b
) else if "%~1"=="run" (
    goto :run
) else if "%~1"=="report" (
    goto :report
) else if "%~1"=="clean" (
    goto :clean
)

echo %BUILD_DIR% created.
@REM Compile each problem into its own executable
:run
for %%p in (%PROBLEMS%) do (
    echo Compiling %%p...
    %CXX% %CXXFLAGS% %SRC_DIR%\%%p.cpp -o %BUILD_DIR%\%%p.exe
    if errorlevel 1 (
        echo Compilation failed for %%p, stopping...
        exit /b 1
    )
)
echo Compilation completed.

for %%p in (%PROBLEMS%) do (
    echo Running %%p.exe...
    call %BUILD_DIR%\%%p.exe
    if errorlevel 1 (
        echo Error running %%p.exe, stopping...
        exit /b 1
    )
)
exit /b

@REM Generate LaTeX report
:report
echo Compiling LaTeX report...
%TEX% %TEX_SRC%
%TEX% %TEX_SRC%  
del /q %TEX_SRC:.tex=.aux%
del /q %TEX_SRC:.tex=.log%


echo Do you want to clean up the intermediate files? [y/n]
set /p clean=
if /i "%clean%"=="y" goto clean
if /i "%clean%"=="n" exit /b 0
exit /b

:clean
echo Cleaning up...
del /q %BUILD_DIR%\*.exe
del /q %TEX_SRC:.tex=.pdf%
echo Clean up completed.
exit /b 0
