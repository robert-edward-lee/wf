#ifndef H_WF
#define H_WF

#include <stddef.h>

#ifndef WF_TYPE
#define WF_TYPE double
#endif /* WF_TYPE */

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
/*                              B-spline windows                              */
/******************************************************************************/
void wf_rect(WF_TYPE *win, size_t N);
void wf_bartlett(WF_TYPE *win, size_t N);
void wf_triang(WF_TYPE *win, size_t N);
void wf_parzen(WF_TYPE *win, size_t N);

/******************************************************************************/
/*                             Cosine-sum windows                             */
/******************************************************************************/
void wf_cosine(WF_TYPE *win, size_t N);
void wf_bohman(WF_TYPE *win, size_t N);
void wf_cosine_sum(WF_TYPE *win, size_t N, const double *a, size_t K);
void wf_general_hamming(WF_TYPE *win, size_t N, double alpha);
void wf_hamming(WF_TYPE *win, size_t N);
void wf_hann(WF_TYPE *win, size_t N);
void wf_blackman_generic(WF_TYPE *win, size_t N, double alpha);
void wf_blackman(WF_TYPE *win, size_t N);
void wf_nuttall(WF_TYPE *win, size_t N);
void wf_blackmanharris(WF_TYPE *win, size_t N);
void wf_flattop(WF_TYPE *win, size_t N);

/******************************************************************************/
/*                             Adjustable windows                             */
/******************************************************************************/
void wf_gaussian(WF_TYPE *win, size_t N, double sigma);
void wf_tukey(WF_TYPE *win, size_t N, double alpha);
void wf_kaiser(WF_TYPE *win, size_t N, double beta);
void wf_kaiser_bessel_derived(WF_TYPE *win, size_t N, double beta);
void wf_poisson(WF_TYPE *win, size_t N, double tau);

/******************************************************************************/
/*                               Hybrid windows                               */
/******************************************************************************/
void wf_barthann(WF_TYPE *win, size_t N);

/******************************************************************************/
/*                                Other windows                               */
/******************************************************************************/
void wf_lanczos(WF_TYPE *win, size_t N);

/******************************************************************************/
/*                                   Utils                                    */
/******************************************************************************/
void wf_print(const WF_TYPE *win, size_t N);

#ifdef __cplusplus
}
#endif

#endif /* H_WF */
