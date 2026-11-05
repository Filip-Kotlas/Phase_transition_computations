@echo off
rem Launcher for PowerShell / cmd: runs cluster.sh with Git's bash (not the WSL bash).
rem   .\slurm\cluster.cmd submit MyRun
"C:\Program Files\Git\bin\bash.exe" "%~dp0cluster.sh" %*
exit /b %ERRORLEVEL%
