@echo off
set SRC_DIR=src
set DATA_DIR=data
set RELEASE_DIR=release
set FIGURE_DIR=figure
set PLOT_DIR=plot

REM 1. 清理之前的生成结果
echo Cleaning old build files...
if exist %RELEASE_DIR%\*.exe del %RELEASE_DIR%\*.exe
if exist *.txt del *.txt
if exist %FIGURE_DIR%\*.png del %FIGURE_DIR%\*.png

REM 2. 编译 src 文件夹中的所有 .cpp 文件，输出到 release 文件夹
echo Compiling C++ files...
for %%f in (%SRC_DIR%\*.cpp) do (
    g++ %%f -o %RELEASE_DIR%\%%~nf.exe
)

REM 3. 运行所有编译好的程序，将结果保存到 src 文件夹中
echo Running compiled programs...
for %%f in (%RELEASE_DIR%\*.exe) do (
    %%f
)

REM 4. 运行 Python 绘图程序，生成图片保存到 figure 文件夹中
echo Running Python plot scripts...
python %PLOT_DIR%\plot_problemB.py
python %PLOT_DIR%\plot_problemC.py

REM 5. 清理 release 文件夹中的可执行文件和 src 文件夹中的输出文本文件
echo Cleaning up...
if exist %RELEASE_DIR%\*.exe del %RELEASE_DIR%\*.exe
if exist %DATA_DIR%\*.txt del %DATA_DIR%\*.txt

echo Process completed. Results saved in %FIGURE_DIR%

