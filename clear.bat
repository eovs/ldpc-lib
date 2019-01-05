@echo off

rmdir /S /Q vs2005\debug
rmdir /S /Q vs2005\release

rmdir /S /Q vs2010\debug
rmdir /S /Q vs2010\release
rmdir /S /Q vs2010\ipch
rmdir /S /Q vs2010\x64


del /A:H /F /Q vs2005\*.suo
del /F /Q vs2005\*.user
del /F /Q vs2005\*.ncb




