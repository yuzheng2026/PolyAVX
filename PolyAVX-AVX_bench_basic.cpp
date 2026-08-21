// PolyAVX-AVX_bench_basic.cpp — PolyAVX 精度验证 + 性能基准一体化测试 (C++98)
//
// 本文件将全功能精度测试与性能基准合并，使用高精度计时器测量每次调用耗时。
// 支持 Windows (QueryPerformanceCounter) 和 Linux/macOS (clock_gettime)。
//
// 编译命令（Linux）：
//   g++ -O3 -march=native -mfma -std=c++98 PolyAVX-AVX_bench_basic.cpp -o PolyAVX-AVX_bench_basic
// 编译命令（Windows MinGW）：
//   g++ -O3 -march=native -mfma -std=c++98 PolyAVX-AVX_bench_basic.cpp -o PolyAVX-AVX_bench_basic.exe
//
// 作者: yuzheng2026 (与 DeepSeek AI 协作开发)
// 许可证: GNU GPLv3 or any later version

#include "poly_avx.hpp"
#include <iostream>
#include <cmath>
#include <cfloat>

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

using namespace poly_avx;

// ==================== 全局数据与辅助函数 ====================

// 测试多项式
PolyD A, B, A0;
const int N = 8;                   // 形式幂级数截断长度
const int BENCH_REPEATS = 10000;   // 性能测试重复次数
// 整数多项式测试数据（用于 PolyI 精确运算）
PolyI IA, IB, IA0;          // 整数多项式，系数在模数范围内
const int N_INT = 8;        // 截断长度
const int MOD_INT = 998244353; // 与 ntt::MOD 一致

// 初始化整数多项式
void setup_int() {
	IA.data.clear();
	IA.data.push_back(1);   // 常数项为 1，便于 log
	IA.data.push_back(2);
	IA.data.push_back(3);   // IA = 1 + 2x + 3x²
	
	IB.data.clear();
	IB.data.push_back(1);
	IB.data.push_back(1);   // IB = 1 + x
	
	IA0 = IA - PolyI(IA[0]); // IA0 = 0 + 2x + 3x²，常数项为 0，便于 exp
}

// 辅助：检查两个整数多项式是否完全相等（误差应为 0）
double coeff_max_error_int(const PolyI& a, const PolyI& b) {
	double err = 0.0;
	int n = a.size() > b.size() ? a.size() : b.size();
	for (int i = 0; i < n; ++i) {
		long long va = (i < a.size()) ? a[i] : 0;
		long long vb = (i < b.size()) ? b[i] : 0;
		double diff = std::abs((double)(va - vb));
		if (diff > err) err = diff;
	}
	return err;
}

// 辅助：检查整数多项式是否与 [1,0,0,...] 相等
double inf_norm_first_one_int(const PolyI& p) {
	double err = 0.0;
	for (int i = 0; i < p.size(); ++i) {
		long long target = (i == 0) ? 1 : 0;
		double diff = std::abs((double)(p[i] - target));
		if (diff > err) err = diff;
	}
	return err;
}
// 高精度计时：返回当前时间的纳秒数
static double now_ns() {
#ifdef _WIN32
	static LARGE_INTEGER freq = {0};
	if (freq.QuadPart == 0) QueryPerformanceFrequency(&freq);
	LARGE_INTEGER counter;
	QueryPerformanceCounter(&counter);
	return double(counter.QuadPart) * 1e9 / double(freq.QuadPart);
#else
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return double(ts.tv_sec) * 1e9 + double(ts.tv_nsec);
#endif
}

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

// 辅助：计算两个多项式逐系数最大绝对误差
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

// 性能测试计时函数：输出每次调用的微秒数
void bench(const char* name, void (*func)(), int repeats = BENCH_REPEATS) {
	double start = now_ns();
	for (int i = 0; i < repeats; ++i) func();
	double end = now_ns();
	double per_call_us = (end - start) / repeats / 1000.0; // 纳秒转微秒
	std::cout << name << ": " << per_call_us << " us/call (" << repeats << " reps)" << std::endl;
}

// 初始化全局多项式
void setup() {
	A.data.clear();
	A.data.push_back(2.0);
	A.data.push_back(3.0);
	A.data.push_back(1.0);   // A = 2 + 3x + x²
	
	B.data.clear();
	B.data.push_back(1.0);
	B.data.push_back(1.0);   // B = 1 + x
	
	A0 = A - PolyD(A[0]);     // A0 = 0 + 3x + x² (常数项为0)
}

// ==================== 性能测试函数 ====================

void bench_add() { PolyD r = A + B; volatile double d = r[0]; (void)d; }
void bench_sub() { PolyD r = A - B; volatile double d = r[0]; (void)d; }
void bench_mul() { PolyD r = A * B; volatile double d = r[0]; (void)d; }
void bench_div() { PolyD r = A / B; volatile double d = r[0]; (void)d; }
void bench_rem() { PolyD r = A % B; volatile double d = r.size(); (void)d; }
void bench_scalar_mul() { PolyD r = A * 2.5; volatile double d = r[0]; (void)d; }
void bench_scalar_div() { PolyD r = B / 2.0; volatile double d = r[0]; (void)d; }

void bench_deriv() { PolyD r = A.deriv(); volatile double d = r[0]; (void)d; }
void bench_integ() { PolyD r = A.integ(); volatile double d = r[0]; (void)d; }

void bench_inv() { PolyD r = A.inv(N); volatile double d = r[0]; (void)d; }
void bench_log() { PolyD r = A.log(N); volatile double d = r[0]; (void)d; }
void bench_exp() { PolyD r = A.exp(N); volatile double d = r[0]; (void)d; }
void bench_sqrt() { PolyD r = A.sqrt(N); volatile double d = r[0]; (void)d; }
void bench_pow_int() { PolyD r = (PolyD(1.0) + B).pow(3, N); volatile double d = r[0]; (void)d; }
void bench_pow_real() { PolyD r = A.pow(0.5, N); volatile double d = r[0]; (void)d; }

void bench_sin() { PolyD r = poly_sin(A0, N); volatile double d = r[0]; (void)d; }
void bench_cos() { PolyD r = poly_cos(A0, N); volatile double d = r[0]; (void)d; }
void bench_tan() { PolyD r = poly_tan(A0, N); volatile double d = r[0]; (void)d; }
void bench_asin() { PolyD r = poly_asin(A0, N); volatile double d = r[0]; (void)d; }
void bench_acos() { PolyD r = poly_acos(A0, N); volatile double d = r[0]; (void)d; }
void bench_atan() { PolyD r = poly_atan(A0, N); volatile double d = r[0]; (void)d; }

void bench_sinh() { PolyD r = poly_sinh(A0, N); volatile double d = r[0]; (void)d; }
void bench_cosh() { PolyD r = poly_cosh(A0, N); volatile double d = r[0]; (void)d; }
void bench_tanh() { PolyD r = poly_tanh(A0, N); volatile double d = r[0]; (void)d; }
void bench_asinh() { PolyD r = poly_asinh(A0, N); volatile double d = r[0]; (void)d; }
void bench_acosh() {
	PolyD C = PolyD(2.0) + A0;    // 常数项 > 1，满足 acosh 定义域
	PolyD r = poly_acosh(C, N);
	volatile double d = r[0]; (void)d;
}
void bench_atanh() { PolyD r = poly_atanh(A0, N); volatile double d = r[0]; (void)d; }

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

void bench_bessel_J1() {
	PolyD r = poly_bessel_J1(8);
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

void bench_erf_series() {
	PolyD r = poly_erf(8);
	volatile double d = r[0]; (void)d;
}

void bench_erfc() {
	PolyD r = poly_erfc(8);
	volatile double d = r[0]; (void)d;
}

void bench_log1p() {
	PolyD A0; A0.data.push_back(0.0); A0.data.push_back(3.0); A0.data.push_back(1.0);
	PolyD r = poly_log1p(A0, 8);
	volatile double d = r[0]; (void)d;
}

// ==================== 精度测试函数 ====================

void run_accuracy_tests() {
	std::cout.precision(15);
	std::cout << "========== 多项式全家桶 全功能精度测试 ==========\n\n";
	
	std::cout << "A = " << A << "\n";
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
	
	PolyD shiftA = poly_shift(A, 2.0, 5);
	std::cout << "A(x+2) = " << shiftA << "\n";
	
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
	
	PolyD Bx; Bx.data.push_back(0.0); Bx.data.push_back(1.0);
	PolyD comp = poly_composite(A, Bx, 5);
	std::cout << "A(B(x)) with B(x)=x = " << comp << "\n";
	err = coeff_max_error(comp.trunc(std::min(5, A.size())), A);
	std::cout << "Composite error (should be 0) = " << err << "\n";
	
	PolyD rev = poly_reversion(Bx, 5);
	std::cout << "reversion of x = " << rev << "   (expected 0 1 0 0 0)\n";
	
	PolyD erfA = poly_erf(A0, 6);
	std::cout << "erf(A0) = " << erfA << "\n";
	std::cout << "erf(0) = " << erfA[0] << "   (should be near 0)\n";
	
	PolyD J0 = poly_bessel_J0(6);
	std::cout << "Bessel J0 (series) = " << J0 << "\n";
	
	PolyD J1 = poly_bessel_J1(12);
	std::cout << "Bessel J1 (series) = " << J1 << "\n";
	
	double x0 = 0.5;
	double val = 0.0, xp = 1.0;
	for (int i = 0; i < J1.size(); ++i) {
		val += J1[i] * xp;
		xp *= x0;
	}
	double j1_ref = 0.2422684576749373;
	double j1_err = std::abs(val - j1_ref);
	std::cout << "| J1(0.5) - reference | = " << j1_err << "   (ideal ~1e-15)\n";
	
	std::cout << "\n========== 精度测试完成 ==========\n\n";
}
void run_accuracy_tests_polyi() {
	std::cout << "\n========== PolyI 精确整数多项式精度测试 ==========\n";
	std::cout.precision(10);
	
	// 乘法验证：IA * IB
	PolyI prod = (IA * IB).trunc(N_INT);
	PolyI prod_manual = PolyI(ntt::convolution(IA.data, IB.data, N_INT));
	double err = coeff_max_error_int(prod, prod_manual);
	std::cout << "| (IA*IB) - NTT卷积 |_inf = " << err << "   (理想 0)\n";
	
	// 求逆验证：IA * IA.inv(n) == 1 (mod x^n)
	PolyI invA = IA.inv(N_INT);
	PolyI check_inv = (IA * invA).trunc(N_INT);
	err = inf_norm_first_one_int(check_inv);
	std::cout << "| IA * IA.inv - 1 |_inf = " << err << "   (理想 0)\n";
	
	// log验证：exp(log(IA)) == IA
	PolyI logA = IA.log(N_INT);
	std::cout << "log(IA) = " << logA << "\n";
	
	PolyI exp_logA = logA.exp(N_INT);
	std::cout << "exp(log(IA)) = " << exp_logA << "\n";
	std::cout << "IA = " << IA << "\n";
	
	err = coeff_max_error_int(exp_logA, IA);
	std::cout << "| exp(log(IA)) - IA |_inf = " << err << "   (理想 0)\n";
	
	// sqrt验证：sqrt(IA)^2 == IA (需要常数项为二次剩余，1是)
	PolyI sqrtA = IA.sqrt(N_INT);
	if (sqrtA.size() == 0) {
		std::cout << "sqrt(IA) 失败：常数项不是二次剩余？\n";
	} else {
		PolyI check_sqrt = (sqrtA * sqrtA).trunc(N_INT);
		err = coeff_max_error_int(check_sqrt, IA);
		std::cout << "| sqrt(IA)^2 - IA |_inf = " << err << "   (理想 0)\n";
	}
	
	// pow验证：(1+x)^3
	PolyI base;
	base.data.push_back(1);   // 常数项 1
	base.data.push_back(1);   // x 项系数 1
	PolyI pow_check = IB.pow(3, N_INT);
	PolyI pow_manual = (IB * IB * IB).trunc(N_INT);
	err = coeff_max_error_int(pow_check, pow_manual);
	std::cout << "| (1+x)^3 - manual |_inf = " << err << "   (理想 0)\n";
	
	// 多模数精确乘法测试
	std::vector<long long> va, vb;
	va.push_back(1000000); va.push_back(2000000); va.push_back(3000000);
	vb.push_back(4000000); vb.push_back(5000000);
	PolyI a(va), b(vb);
	PolyI c = mul_exact(a, b);  // 外部函数
	std::vector<long long> manual = ntt::convolution_multi_mod(va, vb); // 或者朴素
// 比较 c.data 与 manual
	double err_mul = coeff_max_error_int(c, PolyI(manual));
	std::cout << "| mul_exact - manual |_inf = " << err_mul << "   (理想 0)\n";
	
	std::cout << "========== PolyI 精度测试完成 ==========\n\n";
}
// ==================== 性能测试函数 ====================

void run_benchmarks() {
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
	bench("besselJ0 ", bench_bessel_J0);
	bench("besselJ1 ", bench_bessel_J1);
	bench("interpol ", bench_interpolate);
	bench("eval     ", bench_eval);
	bench("erf_series", bench_erf_series);
	bench("erfc     ", bench_erfc);
	bench("log1p    ", bench_log1p);
	
	std::cout << "\nDone.\n";
}

// 在全局区域添加对应的性能测试函数
void bench_polyi_mul() {
	static PolyI bigA, bigB; // 静态变量避免重复构造
	static bool initialized = false;
	if (!initialized) {
		const int N_BIG = 256;
		std::vector<long long> va(N_BIG), vb(N_BIG);
		for (int i = 0; i < N_BIG; ++i) {
			va[i] = (i * 37 + 11) % 1000;
			vb[i] = (i * 91 + 23) % 1000;
		}
		bigA = PolyI(va);
		bigB = PolyI(vb);
		initialized = true;
	}
	PolyI r = bigA * bigB;
	volatile long long d = r[0]; (void)d;
}

void bench_polyi_inv() {
	static PolyI inv_input;
	static bool initialized = false;
	if (!initialized) {
		inv_input.data.clear();
		inv_input.data.push_back(1);
		inv_input.data.push_back(2);
		inv_input.data.push_back(3);
		initialized = true;
	}
	PolyI r = inv_input.inv(64);
	volatile long long d = r[0]; (void)d;
}
void run_benchmarks_polyi() {
	std::cout << "\nPolyI (NTT) Benchmark (repeats=" << BENCH_REPEATS << ")\n";
	
	// 准备两个较大的整数多项式用于性能测试
	const int N_BIG = 256;
	std::vector<long long> va(N_BIG), vb(N_BIG);
	for (int i = 0; i < N_BIG; ++i) {
		va[i] = (i * 37 + 11) % 1000;
		vb[i] = (i * 91 + 23) % 1000;
	}
	PolyI bigA(va), bigB(vb);
	
	// 乘法性能
	bench("PolyI mul (N=256) ", bench_polyi_mul, BENCH_REPEATS);
	// 求逆性能
	bench("PolyI inv (N=64)  ", bench_polyi_inv, BENCH_REPEATS);
}

// ==================== 主函数 ====================

int main() {
	setup();          // 初始化 PolyD 数据
	setup_int();      // 初始化 PolyI 数据
	
	run_accuracy_tests();       // PolyD 精度测试
	run_accuracy_tests_polyi(); // PolyI 精度测试
	
	run_benchmarks();           // PolyD 性能测试
	run_benchmarks_polyi();     // PolyI 性能测试
	
	return 0;
}
