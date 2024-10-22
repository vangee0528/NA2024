@echo off
set SRC_DIR=src
set DATA_DIR=data
set RELEASE_DIR=release
set FIGURE_DIR=figure
set PLOT_DIR=plot
set REPORT_FILE=report.tex

if "%1" == "run" goto run
if "%1" == "clean" goto clean
if "%1" == "report" goto report

echo Invalid command. Use ./make run, ./make clean, or ./make report.
goto end

:run
REM 1. 清理之前的生成结果
echo Cleaning old build files...
if exist %RELEASE_DIR%\*.exe del %RELEASE_DIR%\*.exe
if exist %DATA_DIR%\*.txt del %DATA_DIR%\*.txt
if exist %FIGURE_DIR%\*.png del %FIGURE_DIR%\*.png
if exist *.aux del *.aux
if exist *.log del *.log
if exist *.out del *.out
if exist *.pdf del *.pdf

REM 2. 编译 src 文件夹中的所有 .cpp 文件，输出到 release 文件夹
echo Compiling C++ files...
for %%f in (%SRC_DIR%\*.cpp) do (
    echo Compiling %%f...
    g++ %%f -o %RELEASE_DIR%\%%~nf.exe
)

REM 3. 运行所有编译好的程序，将结果保存到 data 文件夹中
echo Running compiled programs...
for %%f in (%RELEASE_DIR%\*.exe) do (
    %%f
)

REM 4. 运行 Python 绘图程序，生成图片保存到 figure 文件夹中
echo Running Python plot scripts...
python %PLOT_DIR%\plot_problemB.py
python %PLOT_DIR%\plot_problemC.py
python %PLOT_DIR%\plot_problemD.py
python %PLOT_DIR%\plot_problemE.py

echo Process completed. Results saved in %FIGURE_DIR%.
goto end

:clean
REM 清理中间文件
echo Cleaning up...
if exist %RELEASE_DIR%\*.exe del %RELEASE_DIR%\*.exe
if exist %DATA_DIR%\*.txt del %DATA_DIR%\*.txt
if exist %FIGURE_DIR%\*.png del %FIGURE_DIR%\*.png
echo Clean completed.
goto end

:report
REM 编译 LaTeX 文件
echo Compiling LaTeX report...
pdflatex %REPORT_FILE%
echo Report compilation completed.
if exist *.aux del *.aux
if exist *.log del *.log\
if exist *.out del *.out
goto end

:end