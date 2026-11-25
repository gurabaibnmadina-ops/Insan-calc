# Insan Calc

[![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C?logo=cplusplus&logoColor=white)](https://isocpp.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20WSL-lightgrey)](https://github.com/gurabaibnmadina-ops/insan-calc)

Lightweight native **C++17 engineering calculation core** designed for fast workshop and CAD/CAM computations, usable both as a standalone CLI executable and as a shared library (`.so`).

## Project Structure

```text
insan-calc/
├── insan_calc_core.cpp   # C++17 calculation core
├── LICENSE               # MIT License
└── README.md
```

## Build Instructions (Linux / WSL)

### Standalone Executable
```bash
g++ -O3 -march=native -std=c++17 insan_calc_core.cpp -o insan_calc
```

### Shared Library (`.so` for FFI / Python / C Integration)
```bash
g++ -O3 -march=native -std=c++17 -shared -fPIC insan_calc_core.cpp -o libinsan_calc.so
```

## License

Released under the [MIT License](LICENSE).
