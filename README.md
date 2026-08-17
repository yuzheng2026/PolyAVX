# PolyAVX – High‑Performance Polynomial Family Bucket (C++98)

[![License: GPL v3+](https://img.shields.io/badge/License-GPL%20v3%2B-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub release](https://img.shields.io/github/v/release/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX/releases)
[![GitHub issues](https://img.shields.io/github/issues/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX/issues)
[![GitHub pull requests](https://img.shields.io/github/issues-pr/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX/pulls)
[![Last commit](https://img.shields.io/github/last-commit/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX/commits/main)
[![Top language](https://img.shields.io/github/languages/top/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX)
[![Repo size](https://img.shields.io/github/repo-size/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX)
[![CI](https://github.com/yuzheng2026/PolyAVX/actions/workflows/ci.yml/badge.svg)](https://github.com/yuzheng2026/PolyAVX/actions/workflows/ci.yml)
[![PolyChain](https://img.shields.io/badge/Related-PolyChain-blue?style=flat&logo=github)](https://github.com/yuzheng2026/PolyChain)
A single‑header, high‑speed polynomial / formal power series library for **x86‑64** with
**SSE3 / AVX / AVX‑512** and **FMA3** acceleration.
All operations use FFT‑based convolution and are optimised for **double** precision.

A Chinese version of this document is available at [README_CN.md](README_CN.md).

## Features

### Basic algebra
- `+`, `-`, `*` (polynomial and scalar), `/` (polynomial division), `%` (remainder)
- compound assignments `+=`, `-=`, `*=`, `/=`
- equality / inequality operators

### Calculus
- `deriv()` – formal derivative
- `integ()` – formal integral (constant term = 0)

### Formal power series (all modulo `xⁿ`)
- `inv(n)` – multiplicative inverse
- `log(n)` – logarithm (auto‑scaled if constant term ≠ 1)
- `exp(n)` – exponential (auto‑scaled if constant term ≠ 0)
- `sqrt(n)` – square root
- `pow(k, n)` – integer fast‑exponentiation **and** real exponent (via exp‑log)
- `divmod(rhs)` – division with remainder (full polynomials)

### Trigonometric & inverse trigonometric
- `poly_sin(A, n)`, `poly_cos(A, n)`, `poly_tan(A, n)`
- `poly_asin(A, n)`, `poly_acos(A, n)`, `poly_atan(A, n)`

### Hyperbolic & inverse hyperbolic
- `poly_sinh(A, n)`, `poly_cosh(A, n)`, `poly_tanh(A, n)`
- `poly_asinh(A, n)`, `poly_acosh(A, n)`, `poly_atanh(A, n)`
  (integral‑form, log‑identity, and log1p‑stabilised versions available)

### Extended operations
- `poly_shift(A, c, n)` – Taylor shift (`A(x+c)`)
- `multipoint_eval_naive(P, pts)` – evaluate polynomial at multiple points (O(n²))
- `multipoint_interpolate(x, y)` – Lagrange interpolation (O(n²))
- `poly_composite(A, B, n)` – composition `A(B(x))` (Brent‑Kung algorithm)
- `poly_reversion(F, n)` – reversion / compositional inverse (Newton iteration)
- `poly_erf(n)` / `poly_erf(A, n)` – error function (series / composition)
- `poly_bessel_J0(n)` – Bessel J₀ series

### I/O
- `operator<<` and `operator>>` for easy stream input/output of coefficients.

## Accuracy

All tests pass on **GCC 11+** (test environment) and any **C++98 conformant compiler**
including older GCC (3.4+), Clang, and MSVC (2005+).  
Typical infinity‑norm errors for the test polynomial `A = 2 + 3x + x²` (truncated to 8 terms):

| Identity                  | Error       |
|---------------------------|-------------|
| `(1/A) * A`               | ~1.7e-15    |
| `sqrt(A)² - A`            | ~1.3e-15    |
| `exp(log(A)) - A`         | ~1.3e-15    |
| `sin² + cos²`             | ~1.4e-12    |
| `cosh² - sinh²`           | ~4.9e-13    |
| `sin(asin(A₀)) - A₀`     | ~7.0e-11    |
| `sinh(asinh(A₀)) - A₀`   | ~1.3e-10    |
| `tanh(atanh(A₀)) - A₀`   | ~1.2e-09    |

Inverse hyperbolic functions (`asinh`, `acosh`, `atanh`) now use **log1p‑stabilised** formulas
to avoid catastrophic cancellation near the origin and near unity. The residual error in
`tanh(atanh(…))` is dominated by the long chain of FFT, Newton iterations, and truncation
operations – a fundamental limit of double‑precision arithmetic.

## Requirements

- **x86‑64 CPU** with at least SSE3 (virtually all modern CPUs).
- **C++98** compiler (GCC, Clang, MSVC). C++11 or later is also fine.
- **FMA3** support is recommended for maximum accuracy (enabled by `-mfma`).

## Compilation

```bash
g++ -O3 -march=native -mfma -std=c++98 your_program.cpp cpu_dispatch.cpp -o your_program
```

If your CPU doesn’t support FMA, drop `-mfma` – the library falls back to add/sub
SIMD intrinsics automatically.

### Precompiled header (optional)

You can speed up compilation of projects that include `poly_avx.hpp` by generating a
**precompiled header** (`.gch`). Use the **exact same compiler options** as when
building your program:

```bash
g++ -O3 -march=native -mfma -std=c++98 poly_avx.hpp -o poly_avx.hpp.gch
```

Place the resulting `poly_avx.hpp.gch` in the same directory as `poly_avx.hpp`.
GCC will automatically use it if the compile flags match.  
**Note:** Do **not** commit the `.gch` file to version control – it is compiler‑specific
and can be regenerated on demand. Add `*.gch` to your `.gitignore`.

## Runtime CPU dispatch

PolyAVX automatically detects the SIMD capabilities of your CPU at startup and selects
the fastest available implementation of the core complex‑multiplication routine. The
same binary runs optimally on all x86‑64 processors – no recompilation needed.

**Supported paths (in priority order):**
1. **AVX‑512F** (if compiled with `-mavx512f` and supported by the CPU)
2. **AVX + FMA3** (if compiled with `-mavx` and supported by the CPU)
3. **SSE3** (fallback, available on all x86‑64 CPUs)

**How it works:**
- Detection uses GCC/Clang's built‑in `__builtin_cpu_supports` function, called once at startup.
- The global function pointer `pointwise_mul` is set to `pointwise_mul_avx512`,
  `pointwise_mul_avx`, or `pointwise_mul_sse3` accordingly.
- All dispatch logic lives in `cpu_dispatch.cpp` and is initialised before `main()` via a
  global constructor – completely transparent to users.

**Compilation with dispatch:**

```bash
g++ -O3 -march=native -mfma -std=c++98 your_program.cpp cpu_dispatch.cpp -o your_program
```

The `-march=native` flag allows the compiler to generate AVX/AVX‑512 code for your
current CPU; the runtime dispatcher ensures that the resulting binary safely falls
back to SSE3 on older hardware.

If you prefer a **static, single‑file build** (no dispatch), simply omit `cpu_dispatch.cpp`
from the command line – the library then reverts to compile‑time macro selection.
You can also merge `cpu_dispatch.cpp` back into `poly_avx.hpp` to restore a true
single‑header library.

## Quick example

```cpp
#include "poly_avx.hpp"
#include <iostream>
using namespace poly_avx;

int main() {
    // Create polynomial 2 + 3x + x²
    PolyD A;
    A.data.push_back(2.0);
    A.data.push_back(3.0);
    A.data.push_back(1.0);

    // Derivative, integral, log, exp
    std::cout << "A' = " << A.deriv() << "\n";
    std::cout << "∫A = " << A.integ() << "\n";
    std::cout << "log(A) = " << A.log(8) << "\n";
    std::cout << "exp(A) = " << A.exp(8) << "\n";

    // Trigonometric functions (use A0 = A - A[0] for zero constant term)
    PolyD A0 = A - PolyD(A[0]);
    std::cout << "sin(A0) = " << poly_sin(A0, 8) << "\n";

    // Composition
    PolyD B; B.data.push_back(0.0); B.data.push_back(1.0); // B = x
    std::cout << "A(B(x)) = " << poly_composite(A, B, 5) << "\n";

    return 0;
}
```
## Benchmark

A self‑contained benchmark suite (`bench_basic.cpp`) measures the throughput of all
core operations. It executes each function several thousand times and reports the
total wall‑clock time, making it easy to track performance changes across commits.

To build and run the benchmark:

```bash
g++ -O3 -march=native -mfma -std=c++98 bench_basic.cpp cpu_dispatch.cpp -o bench_basic
./bench_basic
```
The output lists 33 operations (algebra, calculus, power series, trigonometric,
hyperbolic, and extended functions) with per‑function timings. The CI pipeline
automatically runs this benchmark on every push and pull request.
## Acknowledgements

This library was created through an extensive collaboration between the author
and **DeepSeek AI**. The AI provided initial code drafts, algorithms explanations,
and debugging assistance; the author performed rigorous testing, optimisation,
and finalisation of every feature.

**Special thanks** to ExplodingKonjac for the original libcp library,
which inspired this project. PolyAVX extends the concept with AVX‑512 support,
additional functions, and a self‑contained C++98 single‑header implementation.
The original library is licensed under GPLv3. PolyAVX is an independent
re‑implementation inspired by it, licensed under GPLv3 to stay true to the same spirit.

## Contributing

Contributions are warmly welcomed! If you're interested in pushing PolyAVX
even closer to the metal, here are some concrete directions:

- **C API** (`extern "C"` wrappers around core functions) – enables linking from
  C, Python, Rust, etc.
- **Hand‑tuned assembly / intrinsics** for FFT butterflies or complex multiply
  (especially for Zen 4, Golden Cove, Xeon series, etc.).
- **Non‑STL memory backend** – replace `std::vector` with a custom allocator so
  the library can be used in kernel / embedded contexts.
- **Runtime CPU dispatch** – detect AVX‑512 / FMA at runtime and select the best
  implementation without recompiling. *(Already implemented — improvements welcome!)*
- **Benchmarks and CI** – a reproducible benchmark suite that tracks performance
  improvements.

If you plan to work on these, please open an issue first so we can discuss the
best approach. The main logic resides in `poly_avx.hpp`; functions like
`pointwise_mul`, `fft`, `convolution`, and `Poly::log` / `Poly::exp` are the
most performance‑critical.

## License

This project is licensed under the **GNU General Public License v3.0 (GPLv3)**
or (at your option) any later version.
See the [LICENSE](LICENSE) file for the full text.

© 2026 yuzheng2026. Licensed under GPLv3.
