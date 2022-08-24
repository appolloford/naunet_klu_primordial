#include <math.h>
/* */
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_sparse.h>  // access to sparse SUNMatrix
/* */
#include "naunet_constants.h"
#include "naunet_macros.h"
#include "naunet_ode.h"
#include "naunet_physics.h"

#define IJth(A, i, j) SM_ELEMENT_D(A, i, j)
#define NVEC_CUDA_CONTENT(x) ((N_VectorContent_Cuda)(x->content))
#define NVEC_CUDA_STREAM(x) (NVEC_CUDA_CONTENT(x)->stream_exec_policy->stream())
#define NVEC_CUDA_BLOCKSIZE(x) \
    (NVEC_CUDA_CONTENT(x)->stream_exec_policy->blockSize())
#define NVEC_CUDA_GRIDSIZE(x, n) \
    (NVEC_CUDA_CONTENT(x)->stream_exec_policy->gridSize(n))

/* */

int Fex(realtype t, N_Vector u, N_Vector udot, void *user_data) {
    /* */
    realtype *y            = N_VGetArrayPointer(u);
    realtype *ydot         = N_VGetArrayPointer(udot);
    NaunetData *u_data     = (NaunetData *)user_data;
    // clang-format off
    realtype nH = u_data->nH;
    realtype Tgas = u_data->Tgas;
    realtype zeta = u_data->zeta;
    realtype Av = u_data->Av;
    realtype omega = u_data->omega;
    realtype mu = u_data->mu;
    realtype gamma = u_data->gamma;
    
    
#if (NHEATPROCS || NCOOLPROCS)
    if (mu < 0) mu = GetMu(y);
    if (gamma < 0) gamma = GetGamma(y);
#endif

    // clang-format on

    realtype k[NREACTIONS] = {0.0};
    EvalRates(k, y, u_data);

#if NHEATPROCS
    realtype kh[NHEATPROCS] = {0.0};
    EvalHeatingRates(kh, y, u_data);
#endif

#if NCOOLPROCS
    realtype kc[NCOOLPROCS] = {0.0};
    EvalCoolingRates(kc, y, u_data);
#endif

    // clang-format off
    ydot[IDX_HeI] = 0.0 - k[3]*y[IDX_HeI]*y[IDX_eM] +
        k[4]*y[IDX_HeII]*y[IDX_eM] + k[5]*y[IDX_HeII]*y[IDX_eM];
    ydot[IDX_HeIII] = 0.0 + k[6]*y[IDX_HeII]*y[IDX_eM] -
        k[7]*y[IDX_HeIII]*y[IDX_eM];
    ydot[IDX_HeII] = 0.0 + k[3]*y[IDX_HeI]*y[IDX_eM] -
        k[4]*y[IDX_HeII]*y[IDX_eM] - k[5]*y[IDX_HeII]*y[IDX_eM] -
        k[6]*y[IDX_HeII]*y[IDX_eM] + k[7]*y[IDX_HeIII]*y[IDX_eM];
    ydot[IDX_H2II] = 0.0 + k[11]*y[IDX_HI]*y[IDX_HII] +
        k[12]*y[IDX_HI]*y[IDX_HII] - k[13]*y[IDX_H2II]*y[IDX_HI] +
        k[14]*y[IDX_H2I]*y[IDX_HII] + k[21]*y[IDX_HM]*y[IDX_HII] -
        k[22]*y[IDX_H2II]*y[IDX_eM] - k[23]*y[IDX_H2II]*y[IDX_eM] -
        k[24]*y[IDX_H2II]*y[IDX_HM];
    ydot[IDX_DII] = 0.0 + k[29]*y[IDX_HII]*y[IDX_DI] -
        k[30]*y[IDX_HI]*y[IDX_DII] - k[31]*y[IDX_H2I]*y[IDX_DII] +
        k[32]*y[IDX_HDI]*y[IDX_HII] - k[37]*y[IDX_DII]*y[IDX_eM];
    ydot[IDX_DI] = 0.0 - k[29]*y[IDX_HII]*y[IDX_DI] +
        k[30]*y[IDX_HI]*y[IDX_DII] - k[33]*y[IDX_H2I]*y[IDX_DI] -
        k[34]*y[IDX_H2I]*y[IDX_DI] + k[35]*y[IDX_HDI]*y[IDX_HI] -
        k[36]*y[IDX_DI]*y[IDX_HM] + k[37]*y[IDX_DII]*y[IDX_eM];
    ydot[IDX_HM] = 0.0 + k[8]*y[IDX_HI]*y[IDX_eM] - k[9]*y[IDX_HM]*y[IDX_HI]
        - k[10]*y[IDX_HM]*y[IDX_HI] - k[17]*y[IDX_HM]*y[IDX_eM] -
        k[18]*y[IDX_HM]*y[IDX_HI] - k[19]*y[IDX_HM]*y[IDX_HI] -
        k[20]*y[IDX_HM]*y[IDX_HII] - k[21]*y[IDX_HM]*y[IDX_HII] -
        k[24]*y[IDX_H2II]*y[IDX_HM] - k[36]*y[IDX_DI]*y[IDX_HM];
    ydot[IDX_HDI] = 0.0 + k[31]*y[IDX_H2I]*y[IDX_DII] -
        k[32]*y[IDX_HDI]*y[IDX_HII] + k[33]*y[IDX_H2I]*y[IDX_DI] +
        k[34]*y[IDX_H2I]*y[IDX_DI] - k[35]*y[IDX_HDI]*y[IDX_HI] +
        k[36]*y[IDX_DI]*y[IDX_HM];
    ydot[IDX_HI] = 0.0 - k[0]*y[IDX_HI]*y[IDX_eM] +
        k[1]*y[IDX_HII]*y[IDX_eM] + k[2]*y[IDX_HII]*y[IDX_eM] -
        k[8]*y[IDX_HI]*y[IDX_eM] - k[9]*y[IDX_HM]*y[IDX_HI] -
        k[10]*y[IDX_HM]*y[IDX_HI] - k[11]*y[IDX_HI]*y[IDX_HII] -
        k[12]*y[IDX_HI]*y[IDX_HII] - k[13]*y[IDX_H2II]*y[IDX_HI] +
        k[14]*y[IDX_H2I]*y[IDX_HII] + k[15]*y[IDX_H2I]*y[IDX_eM] +
        k[15]*y[IDX_H2I]*y[IDX_eM] - k[16]*y[IDX_H2I]*y[IDX_HI] +
        k[16]*y[IDX_H2I]*y[IDX_HI] + k[16]*y[IDX_H2I]*y[IDX_HI] +
        k[16]*y[IDX_H2I]*y[IDX_HI] + k[17]*y[IDX_HM]*y[IDX_eM] -
        k[18]*y[IDX_HM]*y[IDX_HI] + k[18]*y[IDX_HM]*y[IDX_HI] +
        k[18]*y[IDX_HM]*y[IDX_HI] - k[19]*y[IDX_HM]*y[IDX_HI] +
        k[19]*y[IDX_HM]*y[IDX_HI] + k[19]*y[IDX_HM]*y[IDX_HI] +
        k[20]*y[IDX_HM]*y[IDX_HII] + k[20]*y[IDX_HM]*y[IDX_HII] +
        k[22]*y[IDX_H2II]*y[IDX_eM] + k[22]*y[IDX_H2II]*y[IDX_eM] +
        k[23]*y[IDX_H2II]*y[IDX_eM] + k[23]*y[IDX_H2II]*y[IDX_eM] +
        k[24]*y[IDX_H2II]*y[IDX_HM] - k[25]*y[IDX_HI]*y[IDX_HI]*y[IDX_HI] -
        k[25]*y[IDX_HI]*y[IDX_HI]*y[IDX_HI] -
        k[25]*y[IDX_HI]*y[IDX_HI]*y[IDX_HI] +
        k[25]*y[IDX_HI]*y[IDX_HI]*y[IDX_HI] -
        k[26]*y[IDX_HI]*y[IDX_HI]*y[IDX_HI] -
        k[26]*y[IDX_HI]*y[IDX_HI]*y[IDX_HI] -
        k[26]*y[IDX_HI]*y[IDX_HI]*y[IDX_HI] +
        k[26]*y[IDX_HI]*y[IDX_HI]*y[IDX_HI] -
        k[27]*y[IDX_H2I]*y[IDX_HI]*y[IDX_HI] -
        k[27]*y[IDX_H2I]*y[IDX_HI]*y[IDX_HI] -
        k[28]*y[IDX_H2I]*y[IDX_HI]*y[IDX_HI] -
        k[28]*y[IDX_H2I]*y[IDX_HI]*y[IDX_HI] + k[29]*y[IDX_HII]*y[IDX_DI] -
        k[30]*y[IDX_HI]*y[IDX_DII] + k[33]*y[IDX_H2I]*y[IDX_DI] +
        k[34]*y[IDX_H2I]*y[IDX_DI] - k[35]*y[IDX_HDI]*y[IDX_HI];
    ydot[IDX_HII] = 0.0 + k[0]*y[IDX_HI]*y[IDX_eM] -
        k[1]*y[IDX_HII]*y[IDX_eM] - k[2]*y[IDX_HII]*y[IDX_eM] -
        k[11]*y[IDX_HI]*y[IDX_HII] - k[12]*y[IDX_HI]*y[IDX_HII] +
        k[13]*y[IDX_H2II]*y[IDX_HI] - k[14]*y[IDX_H2I]*y[IDX_HII] -
        k[20]*y[IDX_HM]*y[IDX_HII] - k[21]*y[IDX_HM]*y[IDX_HII] -
        k[29]*y[IDX_HII]*y[IDX_DI] + k[30]*y[IDX_HI]*y[IDX_DII] +
        k[31]*y[IDX_H2I]*y[IDX_DII] - k[32]*y[IDX_HDI]*y[IDX_HII];
    ydot[IDX_H2I] = 0.0 + k[9]*y[IDX_HM]*y[IDX_HI] +
        k[10]*y[IDX_HM]*y[IDX_HI] + k[13]*y[IDX_H2II]*y[IDX_HI] -
        k[14]*y[IDX_H2I]*y[IDX_HII] - k[15]*y[IDX_H2I]*y[IDX_eM] -
        k[16]*y[IDX_H2I]*y[IDX_HI] + k[24]*y[IDX_H2II]*y[IDX_HM] +
        k[25]*y[IDX_HI]*y[IDX_HI]*y[IDX_HI] +
        k[26]*y[IDX_HI]*y[IDX_HI]*y[IDX_HI] -
        k[27]*y[IDX_H2I]*y[IDX_HI]*y[IDX_HI] +
        k[27]*y[IDX_H2I]*y[IDX_HI]*y[IDX_HI] +
        k[27]*y[IDX_H2I]*y[IDX_HI]*y[IDX_HI] -
        k[28]*y[IDX_H2I]*y[IDX_HI]*y[IDX_HI] +
        k[28]*y[IDX_H2I]*y[IDX_HI]*y[IDX_HI] +
        k[28]*y[IDX_H2I]*y[IDX_HI]*y[IDX_HI] - k[31]*y[IDX_H2I]*y[IDX_DII] +
        k[32]*y[IDX_HDI]*y[IDX_HII] - k[33]*y[IDX_H2I]*y[IDX_DI] -
        k[34]*y[IDX_H2I]*y[IDX_DI] + k[35]*y[IDX_HDI]*y[IDX_HI];
    ydot[IDX_eM] = 0.0 - k[0]*y[IDX_HI]*y[IDX_eM] + k[0]*y[IDX_HI]*y[IDX_eM]
        + k[0]*y[IDX_HI]*y[IDX_eM] - k[1]*y[IDX_HII]*y[IDX_eM] -
        k[2]*y[IDX_HII]*y[IDX_eM] - k[3]*y[IDX_HeI]*y[IDX_eM] +
        k[3]*y[IDX_HeI]*y[IDX_eM] + k[3]*y[IDX_HeI]*y[IDX_eM] -
        k[4]*y[IDX_HeII]*y[IDX_eM] - k[5]*y[IDX_HeII]*y[IDX_eM] -
        k[6]*y[IDX_HeII]*y[IDX_eM] + k[6]*y[IDX_HeII]*y[IDX_eM] +
        k[6]*y[IDX_HeII]*y[IDX_eM] - k[7]*y[IDX_HeIII]*y[IDX_eM] -
        k[8]*y[IDX_HI]*y[IDX_eM] + k[9]*y[IDX_HM]*y[IDX_HI] +
        k[10]*y[IDX_HM]*y[IDX_HI] - k[15]*y[IDX_H2I]*y[IDX_eM] +
        k[15]*y[IDX_H2I]*y[IDX_eM] - k[17]*y[IDX_HM]*y[IDX_eM] +
        k[17]*y[IDX_HM]*y[IDX_eM] + k[17]*y[IDX_HM]*y[IDX_eM] +
        k[18]*y[IDX_HM]*y[IDX_HI] + k[19]*y[IDX_HM]*y[IDX_HI] +
        k[21]*y[IDX_HM]*y[IDX_HII] - k[22]*y[IDX_H2II]*y[IDX_eM] -
        k[23]*y[IDX_H2II]*y[IDX_eM] + k[36]*y[IDX_DI]*y[IDX_HM] -
        k[37]*y[IDX_DII]*y[IDX_eM];
    ydot[IDX_TGAS] = (gamma - 1.0) * ( 0.0 - kc[0] * y[IDX_HI]*y[IDX_eM] -
        kc[1] * y[IDX_HeI]*y[IDX_eM] - kc[2] * y[IDX_HeII]*y[IDX_eM] - kc[3] *
        y[IDX_HeII]*y[IDX_eM]*y[IDX_eM] - kc[4] * y[IDX_HII]*y[IDX_eM] - kc[5] *
        y[IDX_HeII]*y[IDX_eM] - kc[6] * y[IDX_HeII]*y[IDX_eM] - kc[7] *
        y[IDX_HeIII]*y[IDX_eM] - kc[8] * y[IDX_HI]*y[IDX_eM]*y[IDX_eM] - kc[9] *
        y[IDX_HeII]*y[IDX_eM] - kc[10] * y[IDX_HeII]*y[IDX_eM] ) / kerg /
        GetNumDens(y);
    
    
#if ((NHEATPROCS || NCOOLPROCS) && NAUNET_DEBUG)
    printf("Total heating/cooling rate: %13.7e\n", ydot[IDX_TGAS]);
#endif

    // clang-format on

    /* */

    return NAUNET_SUCCESS;
}