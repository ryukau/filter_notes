Write-Output "`n# cl /nologo /EHsc /std:c++20 /O2"
cl /nologo /EHsc /std:c++20 /O2 /arch:AVX2 benchmark.cpp
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
./benchmark.exe

Write-Output "`n# cl /nologo /EHsc /std:c++20 /O2 /arch:AVX2"
cl /nologo /EHsc /std:c++20 /O2 /arch:AVX2 benchmark.cpp
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
./benchmark.exe

Write-Output "`n# cl /nologo /EHsc /std:c++20 /O2 /arch:AVX2 /fp:fast"
cl /nologo /EHsc /std:c++20 /O2 /arch:AVX2 /fp:fast benchmark.cpp
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
./benchmark.exe

Write-Output "`n"
