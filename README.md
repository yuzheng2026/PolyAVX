# PolyAVX – High‑Performance Polynomial Family Bucket (C++98)

[![License: GPL v3+](https://img.shields.io/badge/License-GPL%20v3%2B-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub release](https://img.shields.io/github/v/release/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX/releases)
[![GitHub issues](https://img.shields.io/github/issues/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX/issues)
[![GitHub pull requests](https://img.shields.io/github/issues-pr/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX/pulls)
[![Last commit](https://img.shields.io/github/last-commit/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX/commits/main)
[![Top language](https://img.shields.io/github/languages/top/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX)
[![Repo size](https://img.shields.io/github/repo-size/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX)

**Related project**: [PolyChain Blockchain Simulator](https://github.com/yuzheng2026/QtWidgetsApplication1)

A high‑speed polynomial / formal power series library for **x86‑64** with **SSE3 / AVX / AVX‑512** and **FMA3** acceleration.  
All operations use FFT‑based convolution and are optimised for **double** precision.  
The library is now **single‑header only** (`poly_avx.hpp`) with built‑in runtime CPU dispatch. No extra `.cpp` files are required.  
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
- `multipoint_eval(P, pts)` – fast multipoint evaluation using a divide‑and‑conquer tree
- `multipoint_interpolate(x, y)` – Lagrange interpolation (O(n²))
- `multipoint_interpolate_fast(x, y)` – fast interpolation (O(n log² n))
- `poly_composite(A, B, n)` – composition `A(B(x))` (Brent‑Kung algorithm)
- `poly_reversion(F, n)` – reversion / compositional inverse (Newton iteration)
- `poly_erf(n)` / `poly_erf(A, n)` – error function (series / composition)
- `poly_erfc(n)` – complementary error function
- `poly_bessel_J0(n)` – Bessel J₀ series
- `poly_bessel_J1(n)` – Bessel J₁ series

### I/O
- `operator<<` and `operator>>` for easy stream input/output of coefficients.

## Accuracy

All identity errors are at or near the limits of double‑precision floating‑point arithmetic.  
Typical test results for the polynomial `A = 2 + 3x + x²` (truncated to 8 terms):

| Identity / Check | Measured Error | Ideal / Expected |
|---|---|---|
| `(1/A)*A - 1` | **0** | ~1e-15 |
| `sqrt(A)² - A` | 4.44089209850063e-16 | ~1e-15 |
| `exp(log(A)) - A` | 1.4658413372004e-16 | ~1e-15 |
| `sin² + cos² - 1` | **0** | ~1e-14 |
| `cosh² - sinh² - 1` | 1.4210854715202e-14 | ~1e-14 |
| `sin(asin(A0)) - A0` | 2.60902410786912e-15 | ~1e-12 |
| `sinh(asinh(A0)) - A0` | 3.37507799486048e-14 | ~1e-13 |
| `tanh(atanh(A0)) - A0` | 1.13686837721616e-13 | ~1e-13 |
| `J1(0.5) - reference` | 6.75015598972095e-14 * | ~1e-15 * |
| Interpolation max error | **0** | <1e-15 |
| Composite error (A(B(x)) with B=x) | **0** | 0 |
| Reversion of x | expected 0 1 0 0 0 | passed |
| `erf(0)` | **0** | near 0 |

\* **Note:** The error of `J1(0.5)` is closely related to the truncation order `n` used in `poly_bessel_J1(n)`.  
For `n=8`, the error is approximately `1.3218556804695e-09`; for `n=12`, it drops to approximately `6.75015598972095e-14`.  
Increase `n` for higher accuracy if needed.

## Performance Comparison

Measured with polynomial truncation length `n=8`, 10000 repetitions per function, using a high-resolution timer.

| Function | Windows (µs) | Linux (µs) |
|----------|--------------|------------|
| `A + B` | 0.0586 | 0.0433 |
| `A - B` | 0.0512 | 0.0250 |
| `A * B` | 0.1578 | 0.1008 |
| `A / B` | 1.1545 | 0.3386 |
| `A % B` | 0.7670 | 0.2315 |
| `A * 2.5` | 0.0501 | 0.0175 |
| `B / 2.0` | 0.0883 | 0.0144 |
| `deriv` | 0.1938 | 0.0303 |
| `integ` | 0.0803 | 0.0197 |
| `inv` | 0.0687 | 0.0337 |
| `log` | 1.2261 | 0.2234 |
| `exp` | 4.2709 | 1.0235 |
| `sqrt` | 0.1456 | 0.0937 |
| `pow(int)` | 1.3473 | 0.3234 |
| `pow(real)` | 4.3581 | 1.5595 |
| `sin` | 5.7077 | 2.0748 |
| `cos` | 5.3239 | 2.1636 |
| `tan` | 5.9694 | 2.2069 |
| `asin` | 7.4112 | 2.9752 |
| `acos` | 7.3243 | 3.0669 |
| `atan` | 0.8993 | 0.2701 |
| `sinh` | 8.8827 | 3.3482 |
| `cosh` | 9.0218 | 3.3697 |
| `tanh` | 9.5994 | 3.5150 |
| `asinh` | 1.0098 | 0.3418 |
| `acosh` | 2.7549 | 0.6498 |
| `atanh` | 1.3221 | 0.2762 |
| `shift` | 0.4815 | 0.1312 |
| `composite` | 2.4084 | 0.5751 |
| `reversion` | 9.0778 | 2.5248 |
| `erf` | 4.0070 | 0.7083 |
| `besselJ0` | 0.1487 | 0.0319 |
| `besselJ1` | 0.1521 | 0.0306 |
| `interpol` | 2.9489 | 0.5897 |
| `eval` | 0.3531 | 0.0968 |
| `erf_series` | 0.1606 | 0.0257 |
| `erfc` | 0.4932 | 0.0649 |
| `log1p` | 1.4894 | 0.3479 |

**Note:** Windows tests used MinGW-w64 GCC 16; Linux tests used WSL2 with GCC 16. Both compiled with `-O3 -march=native -mfma`.

## Requirements

- **x86‑64 CPU** with at least SSE3 (virtually all modern CPUs).
- **C++98** compiler (GCC, Clang, MSVC). C++11 or later is also fine.
- **FMA3** support is recommended for maximum accuracy (enabled by `-mfma`).

## Compilation

Because PolyAVX is now a single‑header library, you only need to include `poly_avx.hpp`:

```bash
g++ -O3 -march=native -mfma -std=c++98 your_program.cpp -o your_program
```

If your CPU doesn’t support FMA, drop `-mfma` – the library falls back to add/sub SIMD intrinsics automatically.  
No separate `.cpp` files or precompiled headers are required.

### Precompiled header (optional)

You can speed up compilation of projects that include `poly_avx.hpp` by generating a **precompiled header** (`.gch`). Use the **exact same compiler options** as when building your program:

```bash
g++ -O3 -march=native -mfma -std=c++98 poly_avx.hpp -o poly_avx.hpp.gch
```

Place the resulting `poly_avx.hpp.gch` in the same directory as `poly_avx.hpp`.  
GCC will automatically use it if the compile flags match.  
**Note:** Do **not** commit the `.gch` file to version control – it is compiler‑specific and can be regenerated on demand. Add `*.gch` to your `.gitignore`.

## Runtime CPU dispatch

PolyAVX automatically detects the SIMD capabilities of your CPU at startup and selects the fastest available implementation of the core complex‑multiplication routine. The same binary runs optimally on all x86‑64 processors – no recompilation needed.

**Supported paths (in priority order):**
1. **AVX‑512F** (if compiled with `-mavx512f` and supported by the CPU)
2. **AVX + FMA3** (if compiled with `-mavx` and supported by the CPU)
3. **SSE3** (fallback, available on all x86‑64 CPUs)

**How it works:**
- Detection uses GCC/Clang's built‑in `__builtin_cpu_supports` function, called once at startup.
- The global function pointer `pointwise_mul` is set to `pointwise_mul_avx512`, `pointwise_mul_avx`, or `pointwise_mul_sse3` accordingly.
- All dispatch logic is now embedded directly inside `poly_avx.hpp`, so no external files are needed.

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

A self‑contained benchmark and accuracy test suite is provided (`PolyAVX-AVX_bench_basic.cpp`). It runs both full functional accuracy tests and performance benchmarks in one program.

To build and run:

```bash
g++ -O3 -march=native -mfma -std=c++98 PolyAVX-AVX_bench_basic.cpp -o PolyAVX-AVX_bench_basic
./PolyAVX-AVX_bench_basic
```

The output lists 33 operations with per‑function timings and verifies key polynomial identities. The CI pipeline can automatically run this program on every push and pull request.

## Acknowledgements

This library was created through an extensive collaboration between the author and **DeepSeek AI**. The AI provided initial code drafts, algorithms explanations, and debugging assistance; the author performed rigorous testing, optimisation, and finalisation of every feature.

**Special thanks** to ExplodingKonjac for the original libcp library, which inspired this project. PolyAVX extends the concept with AVX‑512 support, additional functions, and a self‑contained C++98 single‑header implementation. The original library is licensed under GPLv3. PolyAVX is an independent re‑implementation inspired by it, licensed under GPLv3 to stay true to the same spirit.

## Contributing

Contributions are warmly welcomed! If you're interested in pushing PolyAVX even closer to the metal, here are some concrete directions:

- **C API** (`extern "C"` wrappers around core functions) – enables linking from C, Python, Rust, etc.
- **Hand‑tuned assembly / intrinsics** for FFT butterflies or complex multiply (especially for Zen 4, Golden Cove, Xeon series, etc.).
- **Non‑STL memory backend** – replace `std::vector` with a custom allocator so the library can be used in kernel / embedded contexts.
- **Runtime CPU dispatch** – improve detection logic or add support for more architectures.
- **Benchmarks and CI** – extend the benchmark suite and CI pipeline to cover more edge cases.

If you plan to work on these, please open an issue first so we can discuss the best approach. The main logic resides in `poly_avx.hpp`; functions like `pointwise_mul`, `fft`, `convolution`, and `Poly::log` / `Poly::exp` are the most performance‑critical.

## License

This project is licensed under the **GNU General Public License v3.0 (GPLv3)** or (at your option) any later version.  
See the [LICENSE](LICENSE) file for the full text.

© 2026 yuzheng2026. Licensed under GPLv3.
