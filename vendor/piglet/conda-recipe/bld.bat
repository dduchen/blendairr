:: PIgLET R package build script for Conda (Windows)

"%R%" CMD INSTALL --build .
if errorlevel 1 exit 1
