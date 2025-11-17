# cl /EHsc /W4 /O2 /favor:AMD64 /std:c++20 /fp:precise /DFP_FAST=0 error.cpp
cl /EHsc /W4 /O2 /favor:AMD64 /std:c++20 /fp:fast /DFP_FAST=1 error.cpp
./error.exe
