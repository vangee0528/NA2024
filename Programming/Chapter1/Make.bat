@echo off
rem Define variables
set CXX=g++
set TEX=pdflatex
set SRC_DIR=.
set BUILD_DIR=release
set PROBLEMS=ProblemB ProblemC ProblemD ProblemE
set TEX_SRC=report.tex
set TEX_OUT=report.pdf

rem Create the release directory if it doesn't exist
if not exist %BUILD_DIR% (
    mkdir %BUILD_DIR%
)

rem Compile each problem into its own executable
:compile
for %%p in (%PROBLEMS%) do (
    echo Compiling %%p...
    %CXX% %CXXFLAGS% %SRC_DIR%\%%p.cpp -o %BUILD_DIR%\%%p.exe
    if errorlevel 1 (
        echo Compilation failed for %%p, stopping...
        exit /b 1
    )
)
echo Compilation completed.

rem Run all compiled programs sequentially in the same command prompt
:run
for %%p in (%PROBLEMS%) do (
    echo Running %%p.exe...
    call %BUILD_DIR%\%%p.exe
    if errorlevel 1 (
        echo Error running %%p.exe, stopping...
        exit /b 1
    )
    pause  rem Pause after each program to prevent immediate close
)



rem Compile LaTeX report
:report
echo Compiling LaTeX report...
%TEX% %TEX_SRC%
%TEX% %TEX_SRC%  rem Run twice to ensure cross-references are correct
exit /b 0

rem Check for command line arguments
if "%1" == "run" (
    call :compile
    call :run
) else if "%1" == "clean" (
    call :clean
) else if "%1" == "report" (
    call :report
) else (
    echo Usage: %0 [run|clean|report]
)



exit /b 0
