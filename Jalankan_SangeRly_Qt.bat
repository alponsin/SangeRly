@echo off
rem Menjalankan SangeRly (jendela Qt) dari folder ini. Klik dua kali.
set "SANGERLY_ROOT=%~dp0"
if "%SANGERLY_ROOT:~-1%"=="\" set "SANGERLY_ROOT=%SANGERLY_ROOT:~0,-1%"
cd /d "%~dp0sangerly_qt"
python run.py
if errorlevel 1 pause
