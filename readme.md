# JPEG

A C-based jpeg baseline encoder & decoder.

NTU ITCT HW1

## Prerequirements

- use [bitmap library](https://github.com/ArashPartow/bitmap) to access bmp file.
- this library is written as c++, so we should compile whole project using `g++`

## How to use
- use `makefile` to compile all

```sh
$ make
```

- and you can specific decoder or encoder, even clean (delete)
```sh
$ make encoder
$ make decoder
$ make clean
```

- or compile with `gnu c++`
```sh
$ g++ -o decoder decoder.c -O3
$ g++ -o encoder encoder.c -O3
```

## enviroments

- already tested on my macOS and NTU CSIE server (Arch Linux), so any FreeBSD or Linux  enviroments should be fine.
- more detail below:
```sh
$ g++ -v
Configured with: --prefix=/Applications/Xcode.app/Contents/Developer/usr --with-gxx-include-dir=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk/usr/include/c++/4.2.1
Apple LLVM version 10.0.1 (clang-1001.0.46.4)
Target: x86_64-apple-darwin18.5.0
```

```sh
$ g++ -v

```

## Author

Name: 丁國騰
Student id: R07922009
Department: 資訊工程學系碩士班

## LICENSE

MIT @kuoteng, 2019
