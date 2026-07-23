// cpu_dispatch.cpp
#include "cpu_dispatch.h"
#include <emmintrin.h>
#include <pmmintrin.h>
#ifdef __AVX__
#include <immintrin.h>
#endif

// 跨平台 CPUID
#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <cpuid.h>
#endif

namespace poly_avx {

// ==================== SSE3 版本 ====================
void pointwise_mul_sse3(cd* A, const cd* B, int len) {
    for (int i = 0; i < len; ++i) {
        __m128d a = _mm_loadu_pd((double*)&A[i]);
        __m128d b = _mm_loadu_pd((double*)&B[i]);
#if defined(__FMA__)
        __m128d a_ar = _mm_movedup_pd(a);
        __m128d a_ai = _mm_unpackhi_pd(a, a);
        __m128d b_bi_br = _mm_shuffle_pd(b, b, 1);
        __m128d real_vec = _mm_fnmadd_pd(a_ai, b_bi_br, _mm_mul_pd(a_ar, b));
        __m128d imag_vec = _mm_fmadd_pd(a_ar, b_bi_br, _mm_mul_pd(a_ai, b));
        __m128d res = _mm_shuffle_pd(real_vec, imag_vec, 0);
#else
        __m128d b_sw = _mm_shuffle_pd(b, b, 0x1);
        __m128d m1 = _mm_mul_pd(a, b);
        __m128d m2 = _mm_mul_pd(a, b_sw);
        __m128d X = _mm_unpacklo_pd(m1, m2);
        __m128d Y = _mm_unpackhi_pd(m1, m2);
        __m128d res = _mm_addsub_pd(X, Y);
#endif
        _mm_storeu_pd((double*)&A[i], res);
    }
}

// ==================== AVX 版本 ====================
#ifdef __AVX__
void pointwise_mul_avx(cd* A, const cd* B, int len) {
    for (int i = 0; i <= len - 2; i += 2) {
        __m256d a = _mm256_loadu_pd((double*)&A[i]);
        __m256d b = _mm256_loadu_pd((double*)&B[i]);
#if defined(__FMA__)
        __m256d a_ar = _mm256_permute_pd(a, 0x0);
        __m256d a_ai = _mm256_permute_pd(a, 0xF);
        __m256d b_bi_br = _mm256_permute_pd(b, 0x5);
        __m256d real_vec = _mm256_fnmadd_pd(a_ai, b_bi_br, _mm256_mul_pd(a_ar, b));
        __m256d imag_vec = _mm256_fmadd_pd(a_ar, b_bi_br, _mm256_mul_pd(a_ai, b));
        __m256d res = _mm256_shuffle_pd(real_vec, imag_vec, 0);
#else
        __m256d b_sw = _mm256_permute_pd(b, 0x5);
        __m256d m1 = _mm256_mul_pd(a, b);
        __m256d m2 = _mm256_mul_pd(a, b_sw);
        __m256d X = _mm256_unpacklo_pd(m1, m2);
        __m256d Y = _mm256_unpackhi_pd(m1, m2);
        __m256d res = _mm256_addsub_pd(X, Y);
#endif
        _mm256_storeu_pd((double*)&A[i], res);
    }
    for (int i = len & ~1; i < len; ++i) pointwise_mul_sse3(A + i, B + i, 1);
}
#else
void pointwise_mul_avx(cd*, const cd*, int) { /* 未编译 */ }
#endif

// ==================== AVX-512 版本 ====================
#ifdef __AVX512F__
void pointwise_mul_avx512(cd* A, const cd* B, int len) {
    for (int i = 0; i <= len - 4; i += 4) {
        __m512d a = _mm512_loadu_pd((double*)&A[i]);
        __m512d b = _mm512_loadu_pd((double*)&B[i]);
#if defined(__FMA__)
        __m512d a_ar = _mm512_permute_pd(a, 0x00);
        __m512d a_ai = _mm512_permute_pd(a, 0xFF);
        __m512d b_bi_br = _mm512_permute_pd(b, 0x55);
        __m512d real_vec = _mm512_fnmadd_pd(a_ai, b_bi_br, _mm512_mul_pd(a_ar, b));
        __m512d imag_vec = _mm512_fmadd_pd(a_ar, b_bi_br, _mm512_mul_pd(a_ai, b));
        __m512d res = _mm512_shuffle_pd(real_vec, imag_vec, 0);
#else
        __m512d b_sw = _mm512_permute_pd(b, 0x55);
        __m512d m1 = _mm512_mul_pd(a, b);
        __m512d m2 = _mm512_mul_pd(a, b_sw);
        __m512d X = _mm512_unpacklo_pd(m1, m2);
        __m512d Y = _mm512_unpackhi_pd(m1, m2);
        __m512d res = _mm512_addsub_pd(X, Y);
#endif
        _mm512_storeu_pd((double*)&A[i], res);
    }
    for (int i = len & ~3; i < len; ++i) pointwise_mul_sse3(A + i, B + i, 1);
}
#else
void pointwise_mul_avx512(cd*, const cd*, int) { /* 未编译 */ }
#endif

// ==================== CPU 检测与调度 ====================
pointwise_mul_func g_pointwise_mul = pointwise_mul_sse3;

#if defined(_MSC_VER)
#define CPUID(info, func) __cpuid(info, func)
#else
#define CPUID(info, func) __get_cpuid(func, &info[0], &info[1], &info[2], &info[3])
#endif

void init_cpu_dispatch() {
    int info[4];
    // 检测 AVX-512
#ifdef __AVX512F__
    CPUID(info, 7);
    if (info[1] & (1 << 16)) {  // AVX-512F 位
        g_pointwise_mul = pointwise_mul_avx512;
        return;
    }
#endif
    // 检测 AVX
#ifdef __AVX__
    CPUID(info, 1);
    if (info[2] & (1 << 28)) {  // AVX 位
        g_pointwise_mul = pointwise_mul_avx;
        return;
    }
#endif
    // 默认 SSE3
    g_pointwise_mul = pointwise_mul_sse3;
}

// 自动初始化对象，在 main 之前执行
namespace {
    struct AutoInit {
        AutoInit() { init_cpu_dispatch(); }
    } auto_init;
}

} // namespace poly_avx
