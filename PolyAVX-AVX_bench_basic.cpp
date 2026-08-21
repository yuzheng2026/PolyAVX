// PolyAVX-AVX_bench_basic.cpp — PolyAVX 精度验证 + 性能基准一体化测试
//
// 本文件将全功能精度测试 (AVX.cpp) 与性能基准 (bench_basic.cpp) 合并，
// 方便一次性评估库的正确性与运行效率。
//
// 运行方式：
//   g++ -O3 -march=native -mfma -std=c++98 PolyAVX-AVX_bench_basic.cpp -o PolyAVX-AVX_bench_basic.exe
//   .\PolyAVX-AVX_bench_basic.exe
//
// 作者: yuzheng2026 (与 DeepSeek AI 协作开发)
// 许可证: GNU GPLv3 or any later version

#include "poly_avx.hpp"
#include <iostream>
#include <cmath>
#include <cfloat>
#include <ctime>

using namespace poly_avx;

// ==================== 全局数据与辅助函数 ====================

// 测试多项式
PolyD A, B, A0;
const int N = 8;                   // 形式幂级数截断长度
const int BENCH_REPEATS = 10000;   // 性能测试重复次数

// 辅助：计算多项式与首项为1其余为0的目标向量的最大绝对误差
// 用于验证 inv(A)*A、sin²+cos² 等恒等式
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

// 性能测试计时函数
void bench(const char* name, void (*func)(), int repeats = BENCH_REPEATS) {
	clock_t start = clock();
	for (int i = 0; i < repeats; ++i) func();
	clock_t end = clock();
	double elapsed = double(end - start) / CLOCKS_PER_SEC;
	std::cout << name << ": " << elapsed << " sec (" << repeats << " reps)" << std::endl;
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

// ==================== 主函数 ====================

int main() {
	setup();
	
	// 先运行精度测试
	run_accuracy_tests();
	
	// 再运行性能基准
	run_benchmarks();
	
	return 0;
}
