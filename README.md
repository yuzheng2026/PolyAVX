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

All tests pass on **GCC 11+** (test environment) and any **C++98 conformant compiler**
including older GCC (3.4+), Clang, and MSVC (2005+).  
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

    g++ -O3 -march=native -mfma -std=c++98 your_program.cpp -o your_program

If your CPU doesn’t support FMA, drop `-mfma` – the library falls back to add/sub
SIMD intrinsics automatically.

### Precompiled header (optional)

You can speed up compilation of projects that include `poly_avx.hpp` by generating a
**precompiled header** (`.gch`). Use the **exact same compiler options** as when
building your program:

    g++ -O3 -march=native -mfma -std=c++98 poly_avx.hpp -o poly_avx.hpp.gch

Place the resulting `poly_avx.hpp.gch` in the same directory as `poly_avx.hpp`.
GCC will automatically use it if the compile flags match.  
**Note:** Do **not** commit the `.gch` file to version control – it is compiler‑specific
and can be regenerated on demand. Add `*.gch` to your `.gitignore`.

## Complete example

    #include "poly_avx.hpp"
    #include <iostream>
    #include <vector>
    using namespace poly_avx;

    int main() {
        // ---- Build polynomials ----
        // A = 2 + 3x + x²
        PolyD A;
        A.data.push_back(2.0);
        A.data.push_back(3.0);
        A.data.push_back(1.0);

        // B = 1 + x
        PolyD B;
        B.data.push_back(1.0);
        B.data.push_back(1.0);

        // A0 = A - 2 (constant term becomes 0 for trig/hyperbolic input)
        PolyD A0 = A - PolyD(2.0);

        std::cout.precision(6);
        std::cout << "A = " << A << "\n";
        std::cout << "B = " << B << "\n\n";

        // ---- Basic algebra ----
        std::cout << "A + B = " << (A + B) << "\n";
        std::cout << "A * B = " << (A * B) << "\n";
        std::cout << "A / B = " << (A / B) << "   (quotient)\n";
        std::cout << "A % B = " << (A % B) << "   (remainder)\n\n";

        // ---- Calculus ----
        std::cout << "A' = " << A.deriv() << "\n";
        std::cout << "∫A = " << A.integ() << "\n\n";

        // ---- Formal power series (8 terms) ----
        int N = 8;
        std::cout << "log(A) = " << A.log(N) << "\n";
        std::cout << "exp(A) = " << A.exp(N) << "\n";
        std::cout << "sqrt(A) = " << A.sqrt(N) << "\n";
        std::cout << "1/A = " << A.inv(N) << "\n";
        std::cout << "A^0.5 = " << A.pow(0.5, N) << "\n\n";

        // ---- Trigonometric functions ----
        std::cout << "sin(A0) = " << poly_sin(A0, N) << "\n";
        std::cout << "cos(A0) = " << poly_cos(A0, N) << "\n";
        std::cout << "asin(A0) = " << poly_asin(A0, N) << "\n";
        std::cout << "atan(A0) = " << poly_atan(A0, N) << "\n\n";

        // ---- Hyperbolic functions ----
        std::cout << "sinh(A0) = " << poly_sinh(A0, N) << "\n";
        std::cout << "cosh(A0) = " << poly_cosh(A0, N) << "\n";
        std::cout << "asinh(A0) = " << poly_asinh(A0, N) << "\n";
        std::cout << "atanh(A0) = " << poly_atanh(A0, N) << "\n\n";

        // ---- Extended operations ----
        // Taylor shift: A(x+2)
        std::cout << "A(x+2) = " << poly_shift(A, 2.0, 5) << "\n";

        // Composition A(B(x))
        std::cout << "A(B(x)) = " << poly_composite(A, B, 5) << "\n";

        // Interpolation example
        std::vector<double> xs, ys;
        xs.push_back(0.0); xs.push_back(1.0); xs.push_back(2.0);
        ys.push_back(1.0); ys.push_back(3.0); ys.push_back(7.0);
        PolyD interp = multipoint_interpolate(xs, ys);
        std::cout << "interpolated (points (0,1),(1,3),(2,7)) = " << interp << "\n";

        // Special functions
        std::cout << "erf(A0) = " << poly_erf(A0, 6) << "\n";
        std::cout << "Bessel J0 = " << poly_bessel_J0(6) << "\n";

        return 0;
    }

## Acknowledgements

This library was created through an extensive collaboration between the author
and **DeepSeek AI**. The AI provided initial code drafts, algorithms explanations,
and debugging assistance; the author performed rigorous testing, optimisation,
and finalisation of every feature.

**Special thanks** to ExplodingKonjac for the original 
libcp library,
which inspired this project. PolyAVX extends the concept with AVX-512 support,
additional functions, and a self-contained C++98 single-header implementation.
The original library is licensed under GPLv3. PolyAVX is an independent re‑implementation inspired by it, licensed under LGPLv3 for greater compatibility.

## Contributing

Contributions are warmly welcomed! If you're interested in pushing PolyAVX
even closer to the metal, here are some concrete directions:

- **C API** (`extern "C"` wrappers around core functions) – enables linking from
  C, Python, Rust, etc.
- **Hand‑tuned assembly / intrinsics** for FFT butterflies or complex multiply
  (especially for Zen 4, Golden Cove, or future x86‑64 u‑archs).
- **Non‑STL memory backend** – replace `std::vector` with a custom allocator so
  the library can be used in kernel / embedded contexts.
- **Runtime CPU dispatch** – detect AVX‑512 / FMA at runtime and select the best
  implementation without recompiling.
- **Benchmarks and CI** – a reproducible benchmark suite that tracks performance
  improvements.

If you plan to work on these, please open an issue first so we can discuss the
best approach. The main logic resides in `poly_avx.hpp`; functions like
`pointwise_mul`, `fft`, `convolution`, and `Poly::log` / `Poly::exp` are the
most performance‑critical.

## License
This project is licensed under the **GNU Lesser General Public License v3.0 (LGPLv3)**
or (at your option) any later version.
See the [LICENSE](LICENSE) file for the full text.

© 2026 yuzheng2026. Licensed under LGPLv3.
