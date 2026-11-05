@echo off
rem Same as Sangerly.vbs but keeps the console open so that you can read any error message.
cd /d "%~dp0"
set R_LIBS=
set R_LIBS_USER=
set SANGERLY_DESKTOP=1
"%~dp0R\bin\Rscript.exe" --vanilla launch.R
pause
