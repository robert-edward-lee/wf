#ifndef H_WF
#define H_WF

#define _GNU_SOURCE
#define _USE_MATH_DEFINES

#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif /* M_PI */

#ifndef WF_TYPE
#define WF_TYPE double
#endif /* WF_TYPE */
#ifndef WF_FMT
#define WF_FMT ".05f"
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

/******************************************************************************/
/*                              B-spline windows                              */
/******************************************************************************/
void wf_rect(WF_TYPE *win, size_t N);
void wf_bartlett(WF_TYPE *win, size_t N);
void wf_triang(WF_TYPE *win, size_t N);

/******************************************************************************/
/*                             Cosine-sum windows                             */
/******************************************************************************/
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
void wf_tukey(WF_TYPE *win, size_t N, double alpha);

/******************************************************************************/
/*                                   Utils                                    */
/******************************************************************************/
void wf_print(const WF_TYPE *win, size_t N);

#endif /* H_WF */
