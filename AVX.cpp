// AVX.cpp — PolyAVX 全功能精度验证 (C++98)
//
// 本程序验证 PolyAVX 中所有多项式运算的正确性和精度。
// 通过计算各种数学恒等式（如 inv(A)*A = 1、sin²+cos² = 1 等），
// 检查结果与理论值的最大绝对误差，确保库的可靠性和数值稳定性。
//
// 作者: yuzheng2026 (与 DeepSeek AI 协作开发)
// 许可证: GNU GPLv3 or any later version

#include "poly_avx.hpp"
#include <iostream>
#include <cmath>
#include <cfloat>

using namespace poly_avx;

// 辅助函数：计算多项式与首项为1其余为0的目标向量的最大绝对误差
// 用于验证 inv(A)*A 或 sin²+cos² 等恒等式
double inf_norm_first_one(const PolyD& p) {
    double err = 0.0;
    for (int i = 0; i < p.size(); ++i) {
        double target = (i == 0) ? 1.0 : 0.0;
        double diff = std::abs(p[i] - target);
        if (diff > err) err = diff;
    }
    return err;
}

// 辅助函数：计算两个多项式逐系数最大绝对误差
// 用于比较两个多项式是否接近
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

// 主函数：运行所有精度测试
int main() {
    std::cout.precision(15);
    std::cout << "========== 多项式全家桶 全功能测试 ==========\n\n";

    // 构造 A = 2 + 3x + x^2
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

    // 精度验证
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
    PolyD A0 = A - PolyD(A[0]);   // A0 = 0 + 3x + x^2
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
    err = 0.0;
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

    // 泰勒平移
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

    // 贝塞尔 J1
    PolyD J1 = poly_bessel_J1(12);
    std::cout << "Bessel J1 (series) = " << J1 << "\n";

    // 验证 J1(0.5) 的数值精度
    double x0 = 0.5;
    double val = 0.0, xp = 1.0;
    for (int i = 0; i < J1.size(); ++i) {
        val += J1[i] * xp;
        xp *= x0;
    }
    double j1_ref = 0.2422684576749373;
    double j1_err = std::abs(val - j1_ref);
    std::cout << "| J1(0.5) - reference | = " << j1_err
              << "   (ideal ~1e-15)\n";

    std::cout << "\n========== 所有测试完成 ==========\n";
    return 0;
}
