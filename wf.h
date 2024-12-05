#ifndef H_WF
#define H_WF

#define _GNU_SOURCE
#define _USE_MATH_DEFINES

#include <math.h>
#include <stddef.h>

#ifndef WF_TYPE
#define WF_TYPE double
#endif /* WF_TYPE */
#ifndef WF_FMT
#define WF_FMT ".05lf"
#endif /* WF_FMT */
#ifndef WF_COS
#define WF_COS cos
#endif /* WF_COS */
#ifndef WF_SIN
#define WF_SIN sin
#endif /* WF_SIN */
#ifndef WF_ABS
#define WF_ABS fabs
#endif /* WF_ABS */

void wf_rect(WF_TYPE *win, size_t N);
void wf_triang_generic(WF_TYPE *win, size_t N, size_t L);
void wf_triang(WF_TYPE *win, size_t N);
void wf_bartlett(WF_TYPE *win, size_t N);
void wf_fejer(WF_TYPE *win, size_t N);

/******************************************************************************/
/*                             Cosine-sum windows                             */
/******************************************************************************/
void wf_cosine_sum(WF_TYPE *win, size_t N, const double *a, size_t K);
void wf_hamming_generic(WF_TYPE *win, size_t N, double alpha);
void wf_hamming(WF_TYPE *win, size_t N);
void wf_hann(WF_TYPE *win, size_t N);
void wf_blackman_generic(WF_TYPE *win, size_t N, double alpha);
void wf_blackman(WF_TYPE *win, size_t N);
void wf_nuttal(WF_TYPE *win, size_t N);
void wf_blackman_nuttal(WF_TYPE *win, size_t N);
void wf_blackman_harris(WF_TYPE *win, size_t N);
void wf_flap_top(WF_TYPE *win, size_t N);
/******************************************************************************/
/*                           End cosine-sum windows                           */
/******************************************************************************/

/******************************************************************************/
/*                             Adjustable windows                             */
/******************************************************************************/
void wf_tukey(WF_TYPE *win, size_t N, double alpha);
/******************************************************************************/
/*                           End adjustable windows                           */
/******************************************************************************/

/******************************************************************************/
/*                                   Utils                                    */
/******************************************************************************/
void win_print(const WF_TYPE *win, size_t N);

#endif /* H_WF */
