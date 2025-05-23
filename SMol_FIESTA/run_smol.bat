@echo off
REM Activate the conda environment
CALL %USERPROFILE%\mambaforge\Scripts\activate.bat smol_env REM put the name of your conda environment here (e.g., smol_env)

REM Run the pipeline with your config
python -m SMol_FIESTA.BatchRun_Rebind -c "%USERPROFILE%\Path\To\Your\configfile.toml" REM put the path to your config file here (e.g., %USERPROFILE%\Path\To\Your\configfile.toml)

REM Keep window open to show results or errors
pause
