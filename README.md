# PolyAVX – High‑Performance Polynomial Family Bucket (C++98)

A single‑header, high‑speed polynomial / formal power series library for **x86‑64** with
**SSE3 / AVX / AVX‑512** and **FMA3** acceleration.
All operations use FFT‑based convolution and are optimised for **double** precision.

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
  (both integral‑form and log‑identity versions available)

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

All tests pass on **GCC 11+** with `-O3 -march=native -mfma`.  
Typical infinity‑norm errors for the test polynomial `A = 2 + 3x + x²` (truncated to 8 terms):

| Identity                  | Error       |
|---------------------------|-------------|
| `(1/A) * A`               | ~4e-15      |
| `sqrt(A)² - A`            | ~9e-16      |
| `exp(log(A)) - A`         | ~1.5e-15    |
| `sin² + cos²`             | ~6e-13      |
| `cosh² - sinh²`           | ~9e-13      |
| `sin(asin(A₀)) - A₀`     | ~5e-12      |
| `sinh(asinh(A₀)) - A₀`   | ~5e-12      |
| `tanh(atanh(A₀)) - A₀`   | ~2e-10      |

`atanh` uses the log‑identity form and can be further tuned; the integral‑based
version offers ~2.1e‑10 accuracy.

## Requirements

- **x86‑64 CPU** with at least SSE3 (virtually all modern CPUs).
- **C++98** compiler (GCC, Clang, MSVC). C++11 or later is also fine.
- **FMA3** support is recommended for maximum accuracy (enabled by `-mfma`).

## Compilation

```bash
g++ -O3 -march=native -mfma -std=c++98 your_program.cpp -o your_program
```

If your CPU doesn’t support FMA, drop `-mfma` – the library falls back to add/sub
SIMD intrinsics automatically.

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

## Acknowledgements

This library was created through an extensive collaboration between the author
and **DeepSeek AI**. The AI provided initial code drafts, algorithms explanations,
and debugging assistance; the author performed rigorous testing, optimisation,
and finalisation of every feature.

## License

This project is licensed under the **GNU Lesser General Public License v3.0 (LGPLv3)**.
See the [LICENSE](LICENSE) file for the full text.

---

© 2026 yuzheng2026. Licensed under LGPLv3.
```
