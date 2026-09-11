@echo off
rem ============================================================
rem  CNPS 元素循环分析平台启动器
rem  克隆本仓库后：若自动探测不到 R，请把下面的 RSCRIPT 改成
rem  你机器上的 Rscript.exe 完整路径（任意 R >= 4.2 均可）。
rem ============================================================
set "RSCRIPT=%LOCALAPPDATA%\Programs\R\R-4.4.1\bin\Rscript.exe"

rem --- 依次尝试常见 R 安装位置 ---
if not exist "%RSCRIPT%" if exist "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" set "RSCRIPT=C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
if not exist "%RSCRIPT%" for /d %%D in ("C:\Program Files\R\R-*") do if not exist "%RSCRIPT%" if exist "%%D\bin\Rscript.exe" set "RSCRIPT=%%D\bin\Rscript.exe"
if not exist "%RSCRIPT%" set "RSCRIPT=Rscript.exe"

rem UTF-8 本地编码：保证界面中文正常显示
set LC_CTYPE=.UTF-8

echo Starting CNPS Shiny App from %~dp0 ...
start "CNPS App" /min "%RSCRIPT%" -e "shiny::runApp('%~dp0.', host='127.0.0.1', port=3838, launch.browser=TRUE)"
echo 浏览器将自动打开 http://127.0.0.1:3838 （如未打开请手动访问）
timeout /t 3 >nul
