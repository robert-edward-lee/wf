#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "wf.h"

#ifndef WF_PI
#define WF_PI 3.14159265358979323846264338327950288
#define WF_2PI 6.28318530717958647692528676655900576
#endif /* WF_PI */

#ifndef WF_FMT
#define WF_FMT ".05f"
#endif /* WF_FMT */
#ifndef WF_COS
#define WF_COS cos
#endif /* WF_COS */
#ifndef WF_ACOS
#define WF_ACOS acos
#endif /* WF_ACOS */
#ifndef WF_COSH
#define WF_COSH cosh
#endif /* WF_COSH */
#ifndef WF_ACOSH
#if defined(__STRICT_ANSI__)
extern double acosh(double);
#endif
#define WF_ACOSH acosh
#endif /* WF_ACOSH */
#ifndef WF_SIN
#define WF_SIN sin
#endif /* WF_SIN */
#ifndef WF_ABS
#define WF_ABS fabs
#endif /* WF_ABS */
#ifndef WF_EXP
#define WF_EXP exp
#endif /* WF_EXP */
#ifndef WF_SQRT
#define WF_SQRT sqrt
#endif /* WF_SQRT */
#ifndef WF_POW
#define WF_POW pow
#endif /* WF_POW */

#define WF_COUNTOF(a) (sizeof(a) / sizeof(*(a)))

#define WF_MIN(x, y) ((x) < (y) ? (x) : (y))
#define WF_MAX(x, y) ((x) > (y) ? (x) : (y))

#define WF_SINC(x) ((x != 0.0) ? WF_SIN(WF_PI * (x)) / (WF_PI * (x)) : 1.0)

#ifndef WF_BESSEL_I0
/**
 * @brief Вычисление ряда Чебышёва
 *
 * @note Оригинал взят отсюда: https://www.netlib.org/cephes
 */
static double wf_chbevl(double x, const double *coeff, size_t N) {
    size_t n;
    double b0, b1, b2;

    b0 = *coeff++;
    b1 = 0.0;
    n = N - 1;
    do {
        b2 = b1;
        b1 = b0;
        b0 = x * b1 - b2 + *coeff++;
    } while(--n);

    return 0.5 * (b0 - b2);
}
/**
 * @brief Модифицированная функция Бесселя нулевого порядка
 *
 * @note Оригинал взят отсюда: https://www.netlib.org/cephes
 */
static double wf_i0(double x) {
    static const double A[] = {
        -4.41534164647933937950e-18, 3.33079451882223809783e-17,  -2.43127984654795469359e-16,
        1.71539128555513303061e-15,  -1.16853328779934516808e-14, 7.67618549860493561688e-14,
        -4.85644678311192946090e-13, 2.95505266312963983461e-12,  -1.72682629144155570723e-11,
        9.67580903537323691224e-11,  -5.18979560163526290666e-10, 2.65982372468238665035e-9,
        -1.30002500998624804212e-8,  6.04699502254191894932e-8,   -2.67079385394061173391e-7,
        1.11738753912010371815e-6,   -4.41673835845875056359e-6,  1.64484480707288970893e-5,
        -5.75419501008210370398e-5,  1.88502885095841655729e-4,   -5.76375574538582365885e-4,
        1.63947561694133579842e-3,   -4.32430999505057594430e-3,  1.05464603945949983183e-2,
        -2.37374148058994688156e-2,  4.93052842396707084878e-2,   -9.49010970480476444210e-2,
        1.71620901522208775349e-1,   -3.04682672343198398683e-1,  6.76795274409476084995e-1,
    };
    static const double B[] = {
        -7.23318048787475395456e-18, -4.83050448594418207126e-18, 4.46562142029675999901e-17,
        3.46122286769746109310e-17,  -2.82762398051658348494e-16, -3.42548561967721913462e-16,
        1.77256013305652638360e-15,  3.81168066935262242075e-15,  -9.55484669882830764870e-15,
        -4.15056934728722208663e-14, 1.54008621752140982691e-14,  3.85277838274214270114e-13,
        7.18012445138366623367e-13,  -1.79417853150680611778e-12, -1.32158118404477131188e-11,
        -3.14991652796324136454e-11, 1.18891471078464383424e-11,  4.94060238822496958910e-10,
        3.39623202570838634515e-9,   2.26666899049817806459e-8,   2.04891858946906374183e-7,
        2.89137052083475648297e-6,   6.88975834691682398426e-5,   3.36911647825569408990e-3,
        8.04490411014108831608e-1,
    };

    if(x < 0.0) {
        x = -x;
    }
    if(x <= 8.0) {
        return WF_EXP(x) * wf_chbevl((x / 2.0) - 2.0, A, WF_COUNTOF(A));
    }

    return WF_EXP(x) * wf_chbevl(32.0 / x - 2.0, B, WF_COUNTOF(B)) / WF_SQRT(x);
}
#define WF_BESSEL_I0 wf_i0
#endif /* WF_BESSEL_I0 */

static unsigned bitreverse(unsigned n, unsigned size) {
    unsigned ri = 0;
    while(size != 1) {
        ri *= 2;
        ri |= (n & 1);
        n >>= 1;
        size /= 2;
    }
    return ri;
}
/**
 * @brief In place complex radix-2 FFT.
 */
static void fft_radix2(complex double *z, unsigned size) {
    unsigned i, j;
    unsigned num_subffts, size_subfft;
    complex double *ww;

    /* we start with (size / 2) FFTs of 2 base elements. */
    num_subffts = size / 2;
    size_subfft = 2;

    ww = malloc(size / 2 * sizeof(*ww));
    if(!ww) {
        return;
    }

    for(i = 0; i < size / 2; ++i) {
        ww[i] = cexp(-2.0 * WF_PI * I * i / size);
    }
    /* Permute the input elements (bit-reversal of indices). */
    for(i = 0; i < size; ++i) {
        const unsigned ri = bitreverse(i, size);
        if(i < ri) {
            const complex double temp = z[i];
            z[i] = z[ri];
            z[ri] = temp;
        }
    }
    /* Perform FFTs */
    while(num_subffts != 0) {
        for(i = 0; i < num_subffts; ++i) {
            unsigned subfft_offset = size_subfft * i;

            for(j = 0; j < size_subfft / 2; ++j) {
                unsigned target1 = subfft_offset + j;
                unsigned target2 = subfft_offset + j + size_subfft / 2;

                unsigned left = target1;
                unsigned right = target2;

                const unsigned ww_index = (j * num_subffts);

                const complex double w = ww[ww_index];

                const complex double zleft = z[left];
                const complex double w_zright = w * z[right];

                z[target1] = zleft + w_zright;
                z[target2] = zleft - w_zright;
            }
        }

        num_subffts /= 2;
        size_subfft *= 2;
    }
}

static void czt(complex double *z, unsigned n, complex double *ztrans, unsigned m, complex double w, complex double a) {
    unsigned k;

    complex double *zz, *w2;

    /* Determine next-biggest power-of-two that fits the (n + m - 1) entries we need. */
    unsigned fft_size = 1;
    while(fft_size < n + m - 1) {
        fft_size *= 2;
    }

    zz = malloc(fft_size * sizeof(*zz));
    if(!zz) {
        return;
    }
    w2 = malloc(fft_size * sizeof(*w2));
    if(!w2) {
        free(zz);
        return;
    }

    /* Initialize zz */
    for(k = 0; k < fft_size; ++k) {
        if(k < n) {
            const complex double w1 = cpow(w, 0.5 * k * k) / cpow(a, k);
            zz[k] = w1 * z[k];
        } else {
            zz[k] = 0;
        }
    }
    fft_radix2(zz, fft_size);

    for(k = 0; k < fft_size; ++k) {
        if(k < n + m - 1) {
            const int kshift = k - (n - 1);

            w2[k] = cpow(w, -0.5 * kshift * kshift);
        } else {
            w2[k] = 0;
        }
    }
    fft_radix2(w2, fft_size);

    for(k = 0; k < fft_size; ++k) {
        zz[k] *= w2[k];
    }
    fft_radix2(zz, fft_size);

    /* Make an inverse FFT from the forward FFT.
        - scale all elements by 1 / fft_size;
        - reverse elements 1 .. (fft_size - 1).
    */
    for(k = 0; k < fft_size; ++k) {
        zz[k] /= fft_size;
    }
    for(k = 1; k < fft_size - k; ++k) {
        const unsigned kswap = fft_size - k;

        const complex double temp = zz[k];
        zz[k] = zz[kswap];
        zz[kswap] = temp;
    }

    for(k = 0; k < m; ++k) {
        const complex double w3 = cpow(w, (0.5 * k * k));
        ztrans[k] = w3 * zz[n - 1 + k];
    }
}

static void czt_fft(complex double *z, unsigned size) {
    if(size == 0) {
        return;
    }

    if(__builtin_popcount(size) == 1) {
        fft_radix2(z, size);
    } else {
        const complex double w = cexp(-2.0 * WF_PI * I / size);
        const complex double a = 1;

        czt(z, size, z, size, w, a);
    }
}

/******************************************************************************/
/*                              B-spline windows                              */
/******************************************************************************/
void wf_rect(WF_TYPE *win, size_t N) {
    size_t n;

    if(!win || !N) {
        return;
    }

    for(n = 0; n != N; ++n) {
        win[n] = 1.0;
    }
}

void wf_bartlett(WF_TYPE *win, size_t N) {
    size_t n;

    if(!win || !N) {
        return;
    }

    for(n = 0; n != N; ++n) {
        win[n] = 1.0 - WF_ABS((n - (N - 1.0) / 2.0) / ((N - 1.0) / 2.0));
    }
}

void wf_triang(WF_TYPE *win, size_t N) {
    size_t n;

    if(!win || !N) {
        return;
    }

    if(N % 2) {
        for(n = 0; n != N; ++n) {
            win[n] = 1.0 - WF_ABS((n - (N - 1.0) / 2.0) / ((N + 1.0) / 2.0));
        }
    } else {
        for(n = 0; n != N; ++n) {
            win[n] = 1.0 - WF_ABS((n - (N - 1.0) / 2.0) / (N / 2.0));
        }
    }
}

void wf_parzen(WF_TYPE *win, size_t N) {
    size_t n;
    WF_TYPE x, y;

    if(!win || !N) {
        return;
    }

    for(n = 0; n < N; ++n) {
        x = WF_ABS(2.0 * n - (N - 1)) / N;
        y = 1.0 - x;

        x = 1.0 - 6.0 * x * x + 6.0 * x * x * x;
        y = 2.0 * y * y * y;

        win[n] = WF_MIN(x, y);
    }
}

/******************************************************************************/
/*                             Cosine-sum windows                             */
/******************************************************************************/
void wf_cosine(WF_TYPE *win, size_t N) {
    size_t n;

    if(!win || !N) {
        return;
    }

    for(n = 0; n != N; ++n) {
        win[n] = WF_SIN((WF_PI * (n + 0.5)) / N);
    }
}

void wf_bohman(WF_TYPE *win, size_t N) {
    size_t n;
    double fac;

    if(!win || !N) {
        return;
    }

    for(n = 0; n != N; ++n) {
        fac = fabs(2.0 * n / (N - 1.0) - 1.0);
        win[n] = (1.0 - fac) * WF_COS(WF_PI * fac) + (1.0 / WF_PI) * WF_SIN(WF_PI * fac);
    }
}

void wf_cosine_sum(WF_TYPE *win, size_t N, const double *a, size_t K) {
    int sgn;
    size_t n, k;

    if(!win || !N || !a || !K) {
        return;
    }

    for(n = 0; n != N; ++n) {
        win[n] = a[0];
        sgn = 1;
        for(k = 1; k != K; ++k) {
            sgn *= -1;
            win[n] += sgn * a[k] * WF_COS((WF_2PI * k * n) / (N - 1.0));
        }
    }
}

void wf_general_hamming(WF_TYPE *win, size_t N, double alpha) {
    static double a[2];

    a[0] = alpha;
    a[1] = 1.0 - alpha;
    wf_cosine_sum(win, N, a, sizeof(a) / sizeof(*a));
}

void wf_hamming(WF_TYPE *win, size_t N) {
    wf_general_hamming(win, N, 0.54);
}

void wf_hann(WF_TYPE *win, size_t N) {
    static const double a[] = {0.5, 0.5};
    wf_cosine_sum(win, N, a, sizeof(a) / sizeof(*a));
}

void wf_blackman_generic(WF_TYPE *win, size_t N, double alpha) {
    static double a[3];

    a[0] = (1.0 - alpha) / 2.0;
    a[1] = 1.0 / 2.0;
    a[2] = alpha / 2.0;
    wf_cosine_sum(win, N, a, sizeof(a) / sizeof(*a));
}

void wf_blackman(WF_TYPE *win, size_t N) {
    wf_blackman_generic(win, N, 0.16);
}

void wf_nuttall(WF_TYPE *win, size_t N) {
    static const double a[] = {0.3635819, 0.4891775, 0.1365995, 0.0106411};
    wf_cosine_sum(win, N, a, sizeof(a) / sizeof(*a));
}

void wf_blackmanharris(WF_TYPE *win, size_t N) {
    static const double a[] = {0.35875, 0.48829, 0.14128, 0.01168};
    wf_cosine_sum(win, N, a, sizeof(a) / sizeof(*a));
}

void wf_flattop(WF_TYPE *win, size_t N) {
    static const double a[] = {
        0.21557895,
        0.41663158,
        0.277263158,
        0.083578947,
        0.006947368,
    };
    wf_cosine_sum(win, N, a, sizeof(a) / sizeof(*a));
}

/******************************************************************************/
/*                             Adjustable windows                             */
/******************************************************************************/
void wf_gaussian(WF_TYPE *win, size_t N, double sigma) {
    size_t n;
    double x2;

    if(!win || !N) {
        return;
    }

    for(n = 0; n != N; ++n) {
        x2 = (n - (N - 1.0) / 2.0) * (n - (N - 1.0) / 2.0);
        win[n] = WF_EXP(-0.5 * x2 / (sigma * sigma));
    }
}

void wf_tukey(WF_TYPE *win, size_t N, double alpha) {
    size_t n;

    if(!win || !N) {
        return;
    }

    if(alpha <= 0.0) {
        wf_rect(win, N);
        return;
    } else if(alpha >= 1.0) {
        wf_hann(win, N);
        return;
    }

    for(n = 0; n != (size_t)(alpha * (N - 1.0) / 2.0 + 1.0); ++n) {
        win[n] = (1.0 - WF_COS((WF_2PI * n) / (alpha * (N - 1.0)))) / 2.0;
    }
    for(; n != (size_t)((N - 1.0) / 2.0 + 1.0); ++n) {
        win[n] = 1.0;
    }
    for(; n != N; ++n) {
        win[n] = win[N - n - 1];
    }
}

void wf_kaiser(WF_TYPE *win, size_t N, double beta) {
    size_t n;
    double fac;

    if(!win || !N) {
        return;
    }

    for(n = 0; n != N; ++n) {
        fac = 2.0 * n / (N - 1.0) - 1.0;
        win[n] = WF_BESSEL_I0(beta * WF_SQRT(1.0 - fac * fac)) / WF_BESSEL_I0(beta);
    }
}

void wf_kaiser_bessel_derived(WF_TYPE *win, size_t N, double beta) {
    size_t n;

    if(!win || !N) {
        return;
    }

    wf_kaiser(win, N / 2 + 1, beta);

    for(n = 1; n != N / 2 + 1; ++n) {
        win[n] = win[n - 1] + win[n];
    }
    for(n = 0; n != N / 2 + 1; ++n) {
        win[n] = WF_SQRT(win[n] / win[N / 2]);
    }

    for(n = 0; n != N / 2 + 1; ++n) {
        win[N - 1 - n] = win[n];
    }
}

void wf_chebyshev(double *win, size_t N, double alpha) {
    unsigned i, h;
    double maxw;
    const unsigned order = N - 1;
    const double amplification = WF_POW(10.0, WF_ABS(alpha) / 20.0);
    const double beta = WF_COSH(WF_ACOSH(amplification) / order);

    complex double *p;

    p = malloc(N * sizeof(*p));
    if(!p) {
        return;
    }

    if(N % 2) {
        for(i = 0; i < N; ++i) {
            const double x = beta * WF_COS(WF_PI * i / N);

            if(x > 1.0) {
                p[i] = WF_COSH(order * WF_ACOSH(x));
            } else if(x < -1.0) {
                p[i] = WF_COSH(order * WF_ACOSH(-x));
            } else {
                p[i] = WF_COS(order * WF_ACOS(x));
            }
        }

        czt_fft(p, N);

        /* Example: n = 11
            w[0] w[1] w[2] w[3] w[4] w[5] w[6] w[7] w[8] w[9] w[10]
                                    =
            p[5] p[4] p[3] p[2] p[1] p[0] p[1] p[2] p[3] p[4] p[5]
        */
        h = (N - 1) / 2;
        for(i = 0; i < N; ++i) {
            const unsigned j = (i <= h) ? (h - i) : (i - h);

            win[i] = creal(p[j]);
        }
    } else {
        for(i = 0; i < N; ++i) {
            const double x = beta * WF_COS(WF_PI * i / N);

            const complex double z = cexp(WF_PI * I * i / N);

            if(x > 1) {
                p[i] = z * WF_COSH(order * WF_ACOSH(x));
            } else if(x < -1) {
                p[i] = -z * WF_COSH(order * WF_ACOSH(-x));
            } else {
                p[i] = z * WF_COS(order * WF_ACOS(x));
            }
        }

        czt_fft(p, N);

        /* Example: n = 10
            w[0] w[1] w[2] w[3] w[4] w[5] w[6] w[7] w[8] w[9]
                                    =
            p[5] p[4] p[3] p[2] p[1] p[1] p[2] p[3] p[4] p[5]
        */
        h = N / 2;
        for(i = 0; i < N; ++i) {
            const unsigned j = (i < h) ? (h - i) : (i - h + 1);

            win[i] = creal(p[j]);
        }
    }

    maxw = win[0];
    for(i = 1; i < N; ++i) {
        maxw = WF_MAX(maxw, win[i]);
    }

    for(i = 0; i < N; ++i) {
        win[i] /= maxw;
    }
}

void wf_poisson(WF_TYPE *win, size_t N, double tau) {
    size_t n;

    if(!win || !N) {
        return;
    }

    for(n = 0; n != N; ++n) {
        win[n] = WF_EXP(-fabs(n - (N - 1.0) / 2.0) / tau);
    }
}

/******************************************************************************/
/*                               Hybrid windows                               */
/******************************************************************************/
void wf_barthann(WF_TYPE *win, size_t N) {
    size_t n;
    double fac;

    if(!win || !N) {
        return;
    }

    for(n = 0; n != N; ++n) {
        fac = fabs(n / (N - 1.0) - 0.5);
        win[n] = 0.62 - 0.48 * fac + 0.38 * WF_COS(WF_2PI * fac);
    }
}

/******************************************************************************/
/*                                Other windows                               */
/******************************************************************************/
void wf_lanczos(WF_TYPE *win, size_t N) {
    size_t n;

    if(!win || !N) {
        return;
    }

    for(n = 0; n != N; ++n) {
        win[n] = WF_SINC(2.0 * n / (N - 1.0) - 1.0);
    }
}

/******************************************************************************/
/*                                   Utils                                    */
/******************************************************************************/
void wf_print(const WF_TYPE *win, size_t N) {
    size_t n;

    printf("static const WF_TYPE win[] = {");
    for(n = 0; n != N; ++n) {
        printf("%s%" WF_FMT ",", n % 16 ? " " : "\n    ", win[n]);
    }
    printf("\n};\n");
}
