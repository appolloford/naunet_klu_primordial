#include <math.h>
/* */
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_sparse.h>  // access to sparse SUNMatrix
/* */
/*  */
#include "naunet_ode.h"
/*  */
#include "naunet_constants.h"
#include "naunet_macros.h"
#include "naunet_physics.h"

#define IJth(A, i, j) SM_ELEMENT_D(A, i, j)

/* */

int Fex(realtype t, N_Vector u, N_Vector udot, void *user_data) {
    /* */
    realtype *y            = N_VGetArrayPointer(u);
    realtype *ydot         = N_VGetArrayPointer(udot);
    NaunetData *u_data     = (NaunetData *)user_data;
    // clang-format off
    realtype nH = u_data->nH;
    realtype Tgas = u_data->Tgas;
    
    realtype mu = u_data->mu;
    realtype gamma = u_data->gamma;
    
    if (mu < 0) mu = GetMu(y);
    if (gamma < 0) gamma = GetGamma(y);
    // clang-format on

    realtype k[NREACTIONS] = {0.0};
    EvalRates(k, y, u_data);

    // clang-format off
    ydot[IDX_DI] = 0.0 - k[29]*y[IDX_HII]*y[IDX_DI] +
        k[30]*y[IDX_HI]*y[IDX_DII] - k[33]*y[IDX_H2I]*y[IDX_DI] -
        k[34]*y[IDX_H2I]*y[IDX_DI] + k[35]*y[IDX_HDI]*y[IDX_HI] -
        k[36]*y[IDX_DI]*y[IDX_HM] + k[37]*y[IDX_DII]*y[IDX_eM];
    ydot[IDX_DII] = 0.0 + k[29]*y[IDX_HII]*y[IDX_DI] -
        k[30]*y[IDX_HI]*y[IDX_DII] - k[31]*y[IDX_H2I]*y[IDX_DII] +
        k[32]*y[IDX_HDI]*y[IDX_HII] - k[37]*y[IDX_DII]*y[IDX_eM];
    ydot[IDX_GRAINI] = 0.0 + k[1]*y[IDX_HII]*y[IDX_eM] +
        k[2]*y[IDX_HII]*y[IDX_eM] + k[4]*y[IDX_HeII]*y[IDX_eM] +
        k[5]*y[IDX_HeII]*y[IDX_eM] + k[7]*y[IDX_HeIII]*y[IDX_eM] +
        k[8]*y[IDX_HI]*y[IDX_eM] + k[11]*y[IDX_HI]*y[IDX_HII] +
        k[12]*y[IDX_HI]*y[IDX_HII];
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
    ydot[IDX_HM] = 0.0 + k[8]*y[IDX_HI]*y[IDX_eM] - k[9]*y[IDX_HM]*y[IDX_HI]
        - k[10]*y[IDX_HM]*y[IDX_HI] - k[17]*y[IDX_HM]*y[IDX_eM] -
        k[18]*y[IDX_HM]*y[IDX_HI] - k[19]*y[IDX_HM]*y[IDX_HI] -
        k[20]*y[IDX_HM]*y[IDX_HII] - k[21]*y[IDX_HM]*y[IDX_HII] -
        k[24]*y[IDX_H2II]*y[IDX_HM] - k[36]*y[IDX_DI]*y[IDX_HM];
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
    ydot[IDX_H2II] = 0.0 + k[11]*y[IDX_HI]*y[IDX_HII] +
        k[12]*y[IDX_HI]*y[IDX_HII] - k[13]*y[IDX_H2II]*y[IDX_HI] +
        k[14]*y[IDX_H2I]*y[IDX_HII] + k[21]*y[IDX_HM]*y[IDX_HII] -
        k[22]*y[IDX_H2II]*y[IDX_eM] - k[23]*y[IDX_H2II]*y[IDX_eM] -
        k[24]*y[IDX_H2II]*y[IDX_HM];
    ydot[IDX_HDI] = 0.0 + k[31]*y[IDX_H2I]*y[IDX_DII] -
        k[32]*y[IDX_HDI]*y[IDX_HII] + k[33]*y[IDX_H2I]*y[IDX_DI] +
        k[34]*y[IDX_H2I]*y[IDX_DI] - k[35]*y[IDX_HDI]*y[IDX_HI] +
        k[36]*y[IDX_DI]*y[IDX_HM];
    ydot[IDX_HeI] = 0.0 - k[3]*y[IDX_HeI]*y[IDX_eM] +
        k[4]*y[IDX_HeII]*y[IDX_eM] + k[5]*y[IDX_HeII]*y[IDX_eM];
    ydot[IDX_HeII] = 0.0 + k[3]*y[IDX_HeI]*y[IDX_eM] -
        k[4]*y[IDX_HeII]*y[IDX_eM] - k[5]*y[IDX_HeII]*y[IDX_eM] -
        k[6]*y[IDX_HeII]*y[IDX_eM] + k[7]*y[IDX_HeIII]*y[IDX_eM];
    ydot[IDX_HeIII] = 0.0 + k[6]*y[IDX_HeII]*y[IDX_eM] -
        k[7]*y[IDX_HeIII]*y[IDX_eM];
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
    ydot[IDX_TGAS] = (gamma - 1.0) * ( 0.0 - 1.27e-21 * sqrt(y[IDX_TGAS]) /
        (1.0 + sqrt(y[IDX_TGAS]/1e5)) * exp(-1.578091e5/y[IDX_TGAS]) *
        y[IDX_HI]*y[IDX_eM] - 9.38e-22 * sqrt(y[IDX_TGAS]) / (1.0 +
        sqrt(y[IDX_TGAS]/1e5)) * exp(-2.853354e5/y[IDX_TGAS]) *
        y[IDX_HeI]*y[IDX_eM] - 4.95e-22 * sqrt(y[IDX_TGAS]) / (1.0 +
        sqrt(y[IDX_TGAS]/1e5)) * exp(-6.31515e5/y[IDX_TGAS]) *
        y[IDX_HeII]*y[IDX_eM] - 5.01e-27 * pow(y[IDX_TGAS], -0.1687) / (1.0 +
        sqrt(y[IDX_TGAS]/1e5)) * exp(-5.5338e4/y[IDX_TGAS]) *
        y[IDX_HeII]*y[IDX_eM]*y[IDX_eM] - 8.7e-27 * sqrt(y[IDX_TGAS]) *
        pow(y[IDX_TGAS]/1e3, -0.2) / (1.0+pow(y[IDX_TGAS]/1e6, 0.7)) *
        y[IDX_HII]*y[IDX_eM] - 1.24e-13 * pow(y[IDX_TGAS], -1.5) *
        exp(-4.7e5/y[IDX_TGAS]) * (1.0+0.3*exp(-9.4e4/y[IDX_TGAS])) *
        y[IDX_HeII]*y[IDX_eM] - 1.55e-26 * pow(y[IDX_TGAS], 0.3647) *
        y[IDX_HeII]*y[IDX_eM] - 3.48e-26 * sqrt(y[IDX_TGAS]) *
        pow(y[IDX_TGAS]/1e3, -0.2) / (1.0+pow(y[IDX_TGAS]/1e6, 0.7)) *
        y[IDX_HeIII]*y[IDX_eM] - 9.1e-27 * pow(y[IDX_TGAS], -0.1687) /
        (1.0+sqrt(y[IDX_TGAS]/1e5)) * exp(-1.3179e4/y[IDX_TGAS]) *
        y[IDX_HI]*y[IDX_eM]*y[IDX_eM] - 5.54e-17 * pow(y[IDX_TGAS], -.0397) /
        (1.0+sqrt(y[IDX_TGAS]/1e5)) *exp(-4.73638e5/y[IDX_TGAS]) *
        y[IDX_HeII]*y[IDX_eM] - 5.54e-17 * pow(y[IDX_TGAS], -.0397) /
        (1.0+sqrt(y[IDX_TGAS]/1e5)) *exp(-4.73638e5/y[IDX_TGAS]) *
        y[IDX_HeII]*y[IDX_eM] ) / kerg / GetNumDens(y);
    
    
#if NAUNET_DEBUG
    printf("Total heating/cooling rate: %13.7e\n", ydot[IDX_TGAS]);
#endif
// clang-format on

    /* */

    return NAUNET_SUCCESS;
}