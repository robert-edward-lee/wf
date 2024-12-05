#include <stdio.h>

#include "window_functions.h"

void wf_rect(WF_TYPE *win, size_t N) {
    size_t n;

    if(!win || !N) {
        return;
    }

    for(n = 0; n != N; ++n) {
        win[n] = (WF_TYPE)1.0;
    }
}

void wf_triang_generic(WF_TYPE *win, size_t N, size_t L) {
    size_t n;

    if(!win || !N || !L) {
        return;
    }

    for(n = 0; n != N; ++n) {
        win[n] =
            (WF_TYPE)1.0
            - WF_ABS((n - ((WF_TYPE)N - (WF_TYPE)1.0) / (WF_TYPE)2.0) / (((WF_TYPE)L - (WF_TYPE)1.0) / (WF_TYPE)2.0));
    }
}

void wf_triang(WF_TYPE *win, size_t N) {
    wf_triang_generic(win, N, N);
}

void wf_bartlett(WF_TYPE *win, size_t N) {
    wf_triang_generic(win, N, N + 1);
}

void wf_fejer(WF_TYPE *win, size_t N) {
    wf_triang_generic(win, N, N + 2);
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
        win[n] = (WF_TYPE)a[0];
        sgn = 1;
        for(k = 1; k != K; ++k) {
            sgn *= -1;
            win[n] += sgn * a[k] * WF_COS((2 * M_PI * k * n) / (N - 1));
        }
    }
}

void wf_hamming_generic(WF_TYPE *win, size_t N, double alpha) {
    static double a[2];

    a[0] = alpha;
    a[1] = 1.0 - alpha;
    wf_cosine_sum(win, N, a, sizeof(a) / sizeof(*a));
}

void wf_hamming(WF_TYPE *win, size_t N) {
    static const double a[] = {25.0 / 46.0, 21.0 / 46.0};
    wf_cosine_sum(win, N, a, sizeof(a) / sizeof(*a));
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
    static const double a[] = {3969.0 / 9304.0, 1155.0 / 4652.0, 714.0 / 18608.0};
    wf_cosine_sum(win, N, a, sizeof(a) / sizeof(*a));
}

void wf_nuttal(WF_TYPE *win, size_t N) {
    static const double a[] = {0.355768, 0.487396, 0.144232, 0.012604};
    wf_cosine_sum(win, N, a, sizeof(a) / sizeof(*a));
}

void wf_blackman_nuttal(WF_TYPE *win, size_t N) {
    static const double a[] = {0.3635819, 0.4891775, 0.1365995, 0.0106411};
    wf_cosine_sum(win, N, a, sizeof(a) / sizeof(*a));
}

void wf_blackman_harris(WF_TYPE *win, size_t N) {
    static const double a[] = {0.35875, 0.48829, 0.14128, 0.01168};
    wf_cosine_sum(win, N, a, sizeof(a) / sizeof(*a));
}

void wf_flap_top(WF_TYPE *win, size_t N) {
    static const double a[] = {0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368};
    wf_cosine_sum(win, N, a, sizeof(a) / sizeof(*a));
}
/******************************************************************************/
/*                           End cosine-sum windows                           */
/******************************************************************************/

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

    for(n = 0; n != (size_t)(alpha * (N - 1) / 2.0 + 1); ++n) {
        win[n] = (1 - WF_COS((2.0 * M_PI * n) / (alpha * (N - 1)))) / 2;
    }
    for(; n != (size_t)((N - 1) / 2.0 + 1); ++n) {
        win[n] = 1.0;
    }
    for(; n != N; ++n) {
        win[n] = win[N - n - 1];
    }
}
/******************************************************************************/
/*                           End adjustable windows                           */
/******************************************************************************/

/******************************************************************************/
/*                                   Utils                                    */
/******************************************************************************/
void win_print(const WF_TYPE *win, size_t N) {
    size_t n;

    printf("static const WF_TYPE win[] = {");
    for(n = 0; n != N; ++n) {
        printf("%s%" WF_FMT ",", n % 16 ? " " : "\n    ", win[n]);
    }
    printf("\n};\n");
}
