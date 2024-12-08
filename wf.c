#include <stdio.h>

#include "wf.h"

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

/******************************************************************************/
/*                             Cosine-sum windows                             */
/******************************************************************************/
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
            win[n] += sgn * a[k] * WF_COS((2.0 * M_PI * k * n) / (N - 1.0));
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
        win[n] = (1.0 - WF_COS((2.0 * M_PI * n) / (alpha * (N - 1.0)))) / 2.0;
    }
    for(; n != (size_t)((N - 1.0) / 2.0 + 1.0); ++n) {
        win[n] = 1.0;
    }
    for(; n != N; ++n) {
        win[n] = win[N - n - 1];
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
