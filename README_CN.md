# PolyAVX – 高性能多项式全家桶 (C++98)

[![许可证: GPL v3+](https://img.shields.io/badge/License-GPL%20v3%2B-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub 发布](https://img.shields.io/github/v/release/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX/releases)
[![GitHub Issues](https://img.shields.io/github/issues/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX/issues)
[![GitHub PR](https://img.shields.io/github/issues-pr/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX/pulls)
[![最后提交](https://img.shields.io/github/last-commit/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX/commits/main)
[![主要语言](https://img.shields.io/github/languages/top/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX)
[![仓库大小](https://img.shields.io/github/repo-size/yuzheng2026/PolyAVX)](https://github.com/yuzheng2026/PolyAVX)
[![CI](https://github.com/yuzheng2026/PolyAVX/actions/workflows/ci.yml/badge.svg)](https://github.com/yuzheng2026/PolyAVX/actions/workflows/ci.yml)

**相关项目**: [PolyChain 区块链模拟器](https://github.com/yuzheng2026/PolyChain)

一个用于 **x86‑64** 平台的高性能多项式/形式幂级数单头文件库，全面支持 **SSE3 / AVX / AVX‑512** 和 **FMA3** 指令集加速。  
所有运算均采用基于 FFT 的卷积，并针对 **双精度浮点** 进行深度优化。

本文件的英文版本见 [README.md](README.md)。

## 功能特性

### 基本代数运算
- `+`、`-`、`*`（多项式乘法与标量乘法）、`/`（多项式除法）、`%`（取余）
- 复合赋值运算符 `+=`、`-=`、`*=`、`/=`
- 相等与不等运算符

### 微积分
- `deriv()` —— 形式导数
- `integ()` —— 形式积分（常数项为 0）

### 形式幂级数（所有运算均取模 `xⁿ`）
- `inv(n)` —— 乘法逆
- `log(n)` —— 对数（常数项不为 1 时自动归一化）
- `exp(n)` —— 指数（常数项不为 0 时自动提取因子）
- `sqrt(n)` —— 平方根
- `pow(k, n)` —— 整数快速幂 **及** 实数指数幂（通过 exp‑log 方法）
- `divmod(rhs)` —— 多项式带余除法

### 三角函数与反三角函数
- `poly_sin(A, n)`、`poly_cos(A, n)`、`poly_tan(A, n)`
- `poly_asin(A, n)`、`poly_acos(A, n)`、`poly_atan(A, n)`

### 双曲函数与反双曲函数
- `poly_sinh(A, n)`、`poly_cosh(A, n)`、`poly_tanh(A, n)`
- `poly_asinh(A, n)`、`poly_acosh(A, n)`、`poly_atanh(A, n)`
  （积分形式与对数恒等式形式均提供）

### 扩展运算
- `poly_shift(A, c, n)` —— 泰勒平移（计算 `A(x+c)`）
- `multipoint_eval_naive(P, pts)` —— 多点求值（朴素 O(n²)）
- `multipoint_interpolate(x, y)` —— 拉格朗日插值（O(n²)）
- `poly_composite(A, B, n)` —— 形式幂级数复合 `A(B(x))`（Brent–Kung 算法）
- `poly_reversion(F, n)` —— 复合逆（牛顿迭代法）
- `poly_erf(n)` / `poly_erf(A, n)` —— 误差函数（泰勒级数 / 复合）
- `poly_bessel_J0(n)` —— 第一类零阶贝塞尔函数 J₀

### 输入/输出
- 重载了 `operator<<` 和 `operator>>`，方便系数的流式输入输出。

## 精度说明

所有测试均在 **GCC 11+** 及任何 **符合 C++98 标准的编译器**（包括旧版 GCC 3.4+、Clang、MSVC 2005+）上通过。  
对于测试多项式 `A = 2 + 3x + x²`（截断至 8 项），各项恒等式的无穷范数误差如下：

| 恒等式                     | 误差       |
|----------------------------|-------------|
| `(1/A) * A`                | ~4e-15      |
| `sqrt(A)² - A`             | ~9e-16      |
| `exp(log(A)) - A`          | ~1.5e-15    |
| `sin² + cos²`              | ~6e-13      |
| `cosh² - sinh²`            | ~9e-13      |
| `sin(asin(A₀)) - A₀`      | ~5e-12      |
| `sinh(asinh(A₀)) - A₀`    | ~5e-12      |
| `tanh(atanh(A₀)) - A₀`    | ~2e-10      |

`atanh` 采用对数恒等式实现，可进一步调优；积分版本提供约 2.1e‑10 的精度。

## 编译要求

-   **x86‑64 CPU**，至少需要 SSE3 指令集支持（几乎所有现代 CPU 均满足）。
-   **支持 C++98 的编译器**（GCC、Clang、MSVC）。更高标准（C++11 等）同样兼容。
-   推荐开启 **FMA3** 支持以获得最高精度（通过编译选项 `-mfma` 启用）。

## 编译方法

```bash
g++ -O3 -march=native -mfma -std=c++98 your_program.cpp -o your_program
```
若你的 CPU 不支持 FMA3，去掉 -mfma 即可——库会自动回退到加/减 SIMD 指令。

## 预编译头文件（可选）
你可以通过生成预编译头文件（.gch）来加速包含 poly_avx.hpp 的项目的编译过程。
请使用与编译程序完全相同的选项来生成预编译头：

```bash
g++ -O3 -march=native -mfma -std=c++98 poly_avx.hpp -o poly_avx.hpp.gch
```
将生成的 poly_avx.hpp.gch 放在 poly_avx.hpp 同目录下，GCC 会在编译选项匹配时自动使用它。
注意：请勿将 .gch 文件提交到版本控制——它依赖于编译器环境，可按需重新生成。建议在 .gitignore 中添加 *.gch。

## 快速示例
```cpp
#include "poly_avx.hpp"
#include <iostream>
using namespace poly_avx;

int main() {
    // 构造多项式 A = 2 + 3x + x²
    PolyD A;
    A.data.push_back(2.0);
    A.data.push_back(3.0);
    A.data.push_back(1.0);

    // 导数、积分、对数、指数
    std::cout << "A' = " << A.deriv() << "\n";
    std::cout << "∫A = " << A.integ() << "\n";
    std::cout << "log(A) = " << A.log(8) << "\n";
    std::cout << "exp(A) = " << A.exp(8) << "\n";

    // 三角函数（使用 A0 = A - A[0] 使常数项为 0）
    PolyD A0 = A - PolyD(A[0]);
    std::cout << "sin(A0) = " << poly_sin(A0, 8) << "\n";

    // 复合
    PolyD B; B.data.push_back(0.0); B.data.push_back(1.0); // B = x
    std::cout << "A(B(x)) = " << poly_composite(A, B, 5) << "\n";

    return 0;
}
```
## 致谢
本项目由作者与 DeepSeek AI 深入协作完成。AI 提供了初始代码草稿、算法讲解和调试协助；作者进行了严格的测试、精度优化和所有功能的最终定型。

特别感谢 ExplodingKonjac 的 libcp 项目，PolyAVX 在其启发下诞生，并扩展了 AVX-512 支持、更多函数以及独立的 C++98 单头文件实现。
原项目使用 GPLv3 许可证，PolyAVX 是受其启发的独立重新实现，同样使用 GPLv3 以保持理念一致。

## 贡献指南
热忱欢迎社区贡献！如果你有兴趣将 PolyAVX 推向更底层的极致性能，以下是一些具体方向：

**C API（为核心函数提供 extern "C" 封装）—— 使 C、Python、Rust 等语言可以链接调用。**

**手写汇编 / 内联函数优化 —— 针对 FFT 蝶形运算或复数乘法进行微架构级调优（如 Zen 4、Golden Cove、Xeon 系列等）。**

**非 STL 内存后端 —— 用自定义分配器替换 std::vector，使库可用于内核态或嵌入式环境。**

**运行时 CPU 特性调度 —— 在运行时检测 AVX‑512 / FMA 支持，自动选择最优实现，无需重新编译。**

**性能基准测试与 CI —— 建立可复现的基准测试套件和持续集成，持续追踪性能变化。**

如果你计划在这些方向上贡献代码，请先创建一个 Issue 讨论最佳实现方案。
核心逻辑集中在 poly_avx.hpp 中；函数 pointwise_mul、fft、convolution 以及 Poly::log / Poly::exp 是性能最关键的路径。

## 许可证
本项目采用 GNU General Public License v3.0 (GPLv3) 或（任你选择）任何更高版本进行许可。
完整许可证文本请参阅 LICENSE 文件。

© 2026 yuzheng2026. 基于 GPLv3 许可。
