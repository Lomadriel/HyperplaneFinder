@echo off
set TARGET="Visual Studio 15 2017"
set PLATFORM=v141
set CFGS=Debug;RelWithDebInfo;Release
set OUT=VS2017_%PLATFORM%
mkdir %OUT% 2>NUL
cmd /C "pushd %OUT% & cmake ../.. -G %TARGET% -T %PLATFORM% -DCMAKE_CONFIGURATION_TYPES=%CFGS%"
