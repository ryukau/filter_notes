vcpkg を使ったビルド。

```ps1
mkdir build
cd build
cmake -DCMAKE_TOOLCHAIN_FILE="/src/vcpkg/scripts/buildsystems/vcpkg.cmake" ..
```

64 bit Windows 向けのパッケージのインストール。

```ps1
vcpkg install zlib:x64-windows
vcpkg install zlib openssl --triplet x64-windows
```
