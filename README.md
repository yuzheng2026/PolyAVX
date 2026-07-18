```markdown
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
// test_full.cpp — 多项式全家桶完整测试 (C++98)
#include "poly_avx.hpp"
#include <iostream>
#include <cmath>
#include <cfloat>

using namespace poly_avx;

// 辅助：计算多项式与首项为1其余为0的目标向量的最大绝对误差
double inf_norm_first_one(const PolyD& p) {
    double err = 0.0;
    for (int i = 0; i < p.size(); ++i) {
        double target = (i == 0) ? 1.0 : 0.0;
        double diff = std::abs(p[i] - target);
        if (diff > err) err = diff;
    }
    return err;
}

// 辅助：计算两个多项式逐系数最大绝对误差（忽略长度差异）
double coeff_max_error(const PolyD& a, const PolyD& b) {
    double err = 0.0;
    int n = a.size() > b.size() ? a.size() : b.size();
    for (int i = 0; i < n; ++i) {
        double va = (i < a.size()) ? a[i] : 0.0;
        double vb = (i < b.size()) ? b[i] : 0.0;
        double diff = std::abs(va - vb);
        if (diff > err) err = diff;
    }
    return err;
}

int main() {
    std::cout.precision(15);
    std::cout << "========== 多项式全家桶 全功能测试 ==========\n\n";

    // 构造基础多项式 A = 2 + 3x + x^2
    PolyD A;
    A.data.push_back(2.0);
    A.data.push_back(3.0);
    A.data.push_back(1.0);
    std::cout << "A = " << A << "\n";

    PolyD B;
    B.data.push_back(1.0);
    B.data.push_back(1.0);
    std::cout << "B = " << B << "\n\n";

    // ---- 1. 基本代数 ----
    std::cout << "--- 基本代数 ---\n";
    std::cout << "A + B = " << (A + B) << "\n";
    std::cout << "A - B = " << (A - B) << "\n";
    std::cout << "A * B = " << (A * B) << "\n";
    std::cout << "A / B (商) = " << (A / B) << "\n";
    std::cout << "A % B (余) = " << (A % B) << "\n";
    std::cout << "A * 2.5 = " << (A * 2.5) << "\n";
    std::cout << "B / 2.0 = " << (B / 2.0) << "\n\n";

    // ---- 2. 微积分 ----
    std::cout << "--- 微积分 ---\n";
    PolyD dA = A.deriv();
    std::cout << "A' = " << dA << "\n";
    PolyD iA = A.integ();
    std::cout << "∫A dx (c=0) = " << iA << "\n\n";

    // ---- 3. 形式幂级数 ----
    const int N = 8;
    std::cout << "--- 形式幂级数 (前 " << N << " 项) ---\n";
    PolyD logA = A.log(N);
    std::cout << "log(A) = " << logA << "\n";
    PolyD expA = A.exp(N);
    std::cout << "exp(A) = " << expA << "\n";
    PolyD sqrtA = A.sqrt(N);
    std::cout << "sqrt(A) = " << sqrtA << "\n";
    PolyD invA = A.inv(N);
    std::cout << "1/A = " << invA << "\n";
    PolyD powR = A.pow(0.5, N);
    std::cout << "A^0.5 = " << powR << "\n";
    PolyD powI = (PolyD(1.0) + B).pow(3, N);
    std::cout << "(1+B)^3 = " << powI << "\n\n";

    // 精度验证：常见恒等式
    double err;
    PolyD check;

    check = (invA * A).trunc(N);
    err = inf_norm_first_one(check);
    std::cout << "| (1/A)*A - 1 |_inf = " << err << "   (理想 ~1e-15)\n";

    check = (sqrtA * sqrtA).trunc(N);
    err = coeff_max_error(check, A);
    std::cout << "| sqrt(A)^2 - A |_inf = " << err << "   (理想 ~1e-15)\n";

    check = (logA.exp(N)).trunc(N);
    err = coeff_max_error(check, A);
    std::cout << "| exp(log(A)) - A |_inf = " << err << "   (理想 ~1e-15)\n\n";

    // ---- 4. 三角函数/反三角函数 ----
    std::cout << "--- 三角函数/反三角函数 ---\n";
    PolyD A0 = A - PolyD(A[0]);   // A0 = 0 + 3x + x^2 (常数项为0)
    std::cout << "A0 = A - 2 = " << A0 << "\n";

    PolyD sinA = poly_sin(A0, N);
    PolyD cosA = poly_cos(A0, N);
    std::cout << "sin(A0) = " << sinA << "\n";
    std::cout << "cos(A0) = " << cosA << "\n";

    check = (sinA * sinA + cosA * cosA).trunc(N);
    err = inf_norm_first_one(check);
    std::cout << "| sin²+cos² - 1 |_inf = " << err << "   (理想 ~1e-14)\n";

    PolyD tanA = poly_tan(A0, N);
    std::cout << "tan(A0) = " << tanA << "\n";

    PolyD asinA = poly_asin(A0, N);
    PolyD acosA = poly_acos(A0, N);
    PolyD atanA = poly_atan(A0, N);
    std::cout << "asin(A0) = " << asinA << "\n";
    std::cout << "acos(A0) = " << acosA << "\n";
    std::cout << "atan(A0) = " << atanA << "\n";

    check = (poly_sin(asinA, N) - A0).trunc(N);
    err = inf_norm_first_one(check); // 期望全0，故首项目标为0
    err = 0.0; // 重算
    for (int i = 0; i < check.size(); ++i) {
        double d = std::abs(check[i]);
        if (d > err) err = d;
    }
    std::cout << "| sin(asin(A0)) - A0 |_inf = " << err << "   (理想 ~1e-12)\n\n";

    // ---- 5. 双曲函数/反双曲函数 ----
    std::cout << "--- 双曲函数/反双曲函数 ---\n";
    PolyD sinhA = poly_sinh(A0, N);
    PolyD coshA = poly_cosh(A0, N);
    std::cout << "sinh(A0) = " << sinhA << "\n";
    std::cout << "cosh(A0) = " << coshA << "\n";

    check = (coshA * coshA - sinhA * sinhA).trunc(N);
    err = inf_norm_first_one(check);
    std::cout << "| cosh² - sinh² - 1 |_inf = " << err << "   (理想 ~1e-14)\n";

    PolyD tanhA = poly_tanh(A0, N);
    std::cout << "tanh(A0) = " << tanhA << "\n";

    PolyD asinhA = poly_asinh(A0, N);
    PolyD atanhA = poly_atanh(A0, N);
    std::cout << "asinh(A0) = " << asinhA << "\n";
    std::cout << "atanh(A0) = " << atanhA << "\n";

    check = (poly_sinh(asinhA, N) - A0).trunc(N);
    err = 0.0;
    for (int i = 0; i < check.size(); ++i) {
        double d = std::abs(check[i]);
        if (d > err) err = d;
    }
    std::cout << "| sinh(asinh(A0)) - A0 |_inf = " << err << "   (理想 ~1e-13)\n";

    check = (poly_tanh(atanhA, N) - A0).trunc(N);
    err = 0.0;
    for (int i = 0; i < check.size(); ++i) {
        double d = std::abs(check[i]);
        if (d > err) err = d;
    }
    std::cout << "| tanh(atanh(A0)) - A0 |_inf = " << err << "   (理想 ~1e-13)\n\n";

    // ---- 6. 扩展功能 ----
    std::cout << "--- 扩展功能 ---\n";

    // 平移
    PolyD shiftA = poly_shift(A, 2.0, 5);
    std::cout << "A(x+2) = " << shiftA << "\n";

    // 多点求值与插值
    std::vector<double> xs, ys;
    xs.push_back(0.0); xs.push_back(1.0); xs.push_back(2.0);
    ys.push_back(1.0); ys.push_back(3.0); ys.push_back(7.0);
    PolyD interp = multipoint_interpolate(xs, ys);
    std::cout << "Interpolated (points (0,1),(1,3),(2,7)) = " << interp << "\n";

    std::vector<double> vals = multipoint_eval_naive(interp, xs);
    double interp_err = 0.0;
    for (size_t i = 0; i < xs.size(); ++i)
        interp_err = std::max(interp_err, std::abs(vals[i] - ys[i]));
    std::cout << "Interpolation max error = " << interp_err << "   (理想 < 1e-15)\n";

    // 复合
    PolyD Bx; Bx.data.push_back(0.0); Bx.data.push_back(1.0); // B(x)=x
    PolyD comp = poly_composite(A, Bx, 5);
    std::cout << "A(B(x)) with B(x)=x = " << comp << "\n";
    err = coeff_max_error(comp.trunc(std::min(5, A.size())), A);
    std::cout << "Composite error (should be 0) = " << err << "\n";

    // 复合逆
    PolyD rev = poly_reversion(Bx, 5);
    std::cout << "reversion of x = " << rev << "   (expected 0 1 0 0 0)\n";

    // 特殊函数 erf
    PolyD erfA = poly_erf(A0, 6);
    std::cout << "erf(A0) = " << erfA << "\n";
    std::cout << "erf(0) = " << erfA[0] << "   (should be near 0)\n";

    // 贝塞尔 J0
    PolyD J0 = poly_bessel_J0(6);
    std::cout << "Bessel J0 (series) = " << J0 << "\n";

    std::cout << "\n========== 所有测试完成 ==========\n";
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
