@echo off
rem Define variables
set CXX=g++
set TEX=xelatex
set SRC_DIR=.\src
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
)

:report
echo Compiling LaTeX report...
%TEX% %TEX_SRC%
%TEX% %TEX_SRC%  rem Run twice to ensure cross-references are correct

rem 询问是否清理中间文件
echo Do you want to clean up the intermediate files? [y/n]
set /p clean=
if %clean%==y goto clean
if %clean%==n exit /b 0

:clean
echo Cleaning up...
del /q %BUILD_DIR%\*.exe
@REM del /q %TEX_OUT%
del /q %TEX_SRC:.tex=.aux%
del /q %TEX_SRC:.tex=.log%
del /q %TEX_SRC:.tex=.out%
del /q %TEX_SRC:.tex=.toc%
echo Clean up completed.
exit /b 0



echo Usage: %0 [run|clean|report]
exit /b 0
