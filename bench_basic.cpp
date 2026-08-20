// bench_basic.cpp – PolyAVX 全面性能基准测试 (C++98)
#include "poly_avx.hpp"
#include <iostream>
#include <ctime>
#include <cstdlib>

using namespace poly_avx;

// 辅助：计时并输出
void bench(const char* name, void (*func)(), int repeats = 1000) {
    clock_t start = clock();
    for (int i = 0; i < repeats; ++i) func();
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    std::cout << name << ": " << elapsed << " sec (" << repeats << " reps)" << std::endl;
}

// 全局测试数据
PolyD A, B, A0;
const int N = 8;           // 截断长度
const int BENCH_REPEATS = 5000;

void setup() {
    A.data.clear();
    A.data.push_back(2.0);
    A.data.push_back(3.0);
    A.data.push_back(1.0);   // A = 2 + 3x + x^2

    B.data.clear();
    B.data.push_back(1.0);
    B.data.push_back(1.0);   // B = 1 + x

    A0 = A - PolyD(A[0]);     // A0 = 0 + 3x + x^2 (常数项为0)
}

// ---- 基本代数 ----
void bench_add() { PolyD r = A + B; volatile double d = r[0]; (void)d; }
void bench_sub() { PolyD r = A - B; volatile double d = r[0]; (void)d; }
void bench_mul() { PolyD r = A * B; volatile double d = r[0]; (void)d; }
void bench_div() { PolyD r = A / B; volatile double d = r[0]; (void)d; }
void bench_rem() { PolyD r = A % B; volatile double d = r.size(); (void)d; }
void bench_scalar_mul() { PolyD r = A * 2.5; volatile double d = r[0]; (void)d; }
void bench_scalar_div() { PolyD r = B / 2.0; volatile double d = r[0]; (void)d; }

// ---- 微积分 ----
void bench_deriv() { PolyD r = A.deriv(); volatile double d = r[0]; (void)d; }
void bench_integ() { PolyD r = A.integ(); volatile double d = r[0]; (void)d; }

// ---- 形式幂级数 ----
void bench_inv() { PolyD r = A.inv(N); volatile double d = r[0]; (void)d; }
void bench_log() { PolyD r = A.log(N); volatile double d = r[0]; (void)d; }
void bench_exp() { PolyD r = A.exp(N); volatile double d = r[0]; (void)d; }
void bench_sqrt() { PolyD r = A.sqrt(N); volatile double d = r[0]; (void)d; }
void bench_pow_int() { PolyD r = (PolyD(1.0) + B).pow(3, N); volatile double d = r[0]; (void)d; }
void bench_pow_real() { PolyD r = A.pow(0.5, N); volatile double d = r[0]; (void)d; }

// ---- 三角函数 ----
void bench_sin() { PolyD r = poly_sin(A0, N); volatile double d = r[0]; (void)d; }
void bench_cos() { PolyD r = poly_cos(A0, N); volatile double d = r[0]; (void)d; }
void bench_tan() { PolyD r = poly_tan(A0, N); volatile double d = r[0]; (void)d; }
void bench_asin() { PolyD r = poly_asin(A0, N); volatile double d = r[0]; (void)d; }
void bench_acos() { PolyD r = poly_acos(A0, N); volatile double d = r[0]; (void)d; }
void bench_atan() { PolyD r = poly_atan(A0, N); volatile double d = r[0]; (void)d; }

// ---- 双曲函数 ----
void bench_sinh() { PolyD r = poly_sinh(A0, N); volatile double d = r[0]; (void)d; }
void bench_cosh() { PolyD r = poly_cosh(A0, N); volatile double d = r[0]; (void)d; }
void bench_tanh() { PolyD r = poly_tanh(A0, N); volatile double d = r[0]; (void)d; }
void bench_asinh() { PolyD r = poly_asinh(A0, N); volatile double d = r[0]; (void)d; }
void bench_acosh() {
    PolyD C = PolyD(2.0) + A0;    // 常数项 > 1
    PolyD r = poly_acosh(C, N);
    volatile double d = r[0]; (void)d;
}
void bench_atanh() { PolyD r = poly_atanh(A0, N); volatile double d = r[0]; (void)d; }

// ---- 扩展运算 ----
void bench_shift() { PolyD r = poly_shift(A, 2.0, 5); volatile double d = r[0]; (void)d; }

void bench_composite() {
    PolyD Bx; Bx.data.push_back(0.0); Bx.data.push_back(1.0); // Bx = x
    PolyD r = poly_composite(A, Bx, 5);
    volatile double d = r[0]; (void)d;
}

void bench_reversion() {
    PolyD F; F.data.push_back(0.0); F.data.push_back(1.0); // F = x
    PolyD r = poly_reversion(F, 5);
    volatile double d = r[0]; (void)d;
}

void bench_erf() {
    PolyD r = poly_erf(A0, 6);
    volatile double d = r[0]; (void)d;
}

void bench_bessel_J0() {
    PolyD r = poly_bessel_J0(6);
    volatile double d = r[0]; (void)d;
}

void bench_interpolate() {
    std::vector<double> xs, ys;
    xs.push_back(0.0); xs.push_back(1.0); xs.push_back(2.0);
    ys.push_back(1.0); ys.push_back(3.0); ys.push_back(7.0);
    PolyD r = multipoint_interpolate(xs, ys);
    volatile double d = r[0]; (void)d;
}

void bench_eval() {
    PolyD P; P.data.push_back(1.0); P.data.push_back(1.0); P.data.push_back(1.0);
    std::vector<double> pts; pts.push_back(0.0); pts.push_back(1.0); pts.push_back(2.0);
    std::vector<double> vals = multipoint_eval_naive(P, pts);
    volatile double d = vals[0]; (void)d;
}

// ---- 补充：特殊函数 ----
void bench_erf_series() {
    PolyD r = poly_erf(8);               // 直接生成 erf 级数
    volatile double d = r[0]; (void)d;
}

void bench_erfc() {
    PolyD r = poly_erfc(8);              // 生成 erfc 级数
    volatile double d = r[0]; (void)d;
}

void bench_bessel_J1() {
    PolyD r = poly_bessel_J1(8);         // 生成 J1 级数
    volatile double d = r[0]; (void)d;
}

void bench_log1p() {
    PolyD A0;                            // 常数项为 0
    A0.data.push_back(0.0);
    A0.data.push_back(3.0);
    A0.data.push_back(1.0);
    PolyD r = poly_log1p(A0, 8);         // log(1 + A0)
    volatile double d = r[0]; (void)d;
}

int main() {
    setup();
    std::cout << "PolyAVX Benchmark (repeats=" << BENCH_REPEATS << ")\n\n";

    bench("A+B      ", bench_add);
    bench("A-B      ", bench_sub);
    bench("A*B      ", bench_mul);
    bench("A/B      ", bench_div);
    bench("A%B      ", bench_rem);
    bench("A*2.5    ", bench_scalar_mul);
    bench("B/2.0    ", bench_scalar_div);
    bench("deriv    ", bench_deriv);
    bench("integ    ", bench_integ);
    bench("inv      ", bench_inv);
    bench("log      ", bench_log);
    bench("exp      ", bench_exp);
    bench("sqrt     ", bench_sqrt);
    bench("pow(int) ", bench_pow_int);
    bench("pow(real)", bench_pow_real);
    bench("sin      ", bench_sin);
    bench("cos      ", bench_cos);
    bench("tan      ", bench_tan);
    bench("asin     ", bench_asin);
    bench("acos     ", bench_acos);
    bench("atan     ", bench_atan);
    bench("sinh     ", bench_sinh);
    bench("cosh     ", bench_cosh);
    bench("tanh     ", bench_tanh);
    bench("asinh    ", bench_asinh);
    bench("acosh    ", bench_acosh);
    bench("atanh    ", bench_atanh);
    bench("shift    ", bench_shift);
    bench("composite", bench_composite);
    bench("reversion", bench_reversion);
    bench("erf      ", bench_erf);
    bench("besselJ0 ", bench_bessel_J0);        // 原有 J0
    bench("besselJ1 ", bench_bessel_J1);     // 新增 J1
    bench("interpol ", bench_interpolate);
    bench("eval     ", bench_eval);
    bench("erf_series", bench_erf_series);
    bench("erfc     ", bench_erfc);
    bench("log1p    ", bench_log1p);
    
    std::cout << "\nDone.\n";
    return 0;
}
