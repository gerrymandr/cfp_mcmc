rem rd /s /q build
rem mkdir build
rem pushd build

rem set

rem set CFLAGS=
rem set CXXFLAGS=

rem cmake -G"%CMAKE_GENERATOR%" ../
rem if errorlevel 1 exit 1

rem cmake --build . --config Release
rem if errorlevel 1 exit 1

COPY chain.exe %LIBRARY_BIN%

popd
COPY CurrentRep.txt %PREFIX%
xcopy %RECIPE_DIR%\..\bill_plans %PREFIX%\bill_plans /E /i
COPY *.csv %PREFIX%
xcopy %RECIPE_DIR%\..\scripts %PREFIX%\scripts /E /i