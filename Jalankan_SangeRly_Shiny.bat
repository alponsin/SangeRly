@echo off
rem Menjalankan SangeRly versi Shiny (di browser). Klik dua kali.
cd /d "%~dp0"
Rscript -e "shiny::runApp('.', launch.browser = TRUE)"
if errorlevel 1 pause
