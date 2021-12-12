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

int Jac(realtype t, N_Vector u, N_Vector fu, SUNMatrix jmatrix, void *user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    /* */
    realtype *y            = N_VGetArrayPointer(u);
    sunindextype *rowptrs  = SUNSparseMatrix_IndexPointers(jmatrix);
    sunindextype *colvals  = SUNSparseMatrix_IndexValues(jmatrix);
    realtype *data         = SUNSparseMatrix_Data(jmatrix);
    NaunetData *u_data     = (NaunetData *)user_data;

    // clang-format off
    realtype mu = u_data->mu;
    realtype gamma = u_data->gamma;
        
    if (mu < 0) mu = GetMu(y);
    if (gamma < 0) gamma = GetGamma(y);
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
    // number of non-zero elements in each row
    rowptrs[0] = 0;
    rowptrs[1] = 8;
    rowptrs[2] = 15;
    rowptrs[3] = 20;
    rowptrs[4] = 29;
    rowptrs[5] = 38;
    rowptrs[6] = 44;
    rowptrs[7] = 53;
    rowptrs[8] = 59;
    rowptrs[9] = 66;
    rowptrs[10] = 69;
    rowptrs[11] = 73;
    rowptrs[12] = 76;
    rowptrs[13] = 87;
    rowptrs[14] = 93;
    
    // the column index of non-zero elements
    colvals[0] = 0;
    colvals[1] = 1;
    colvals[2] = 3;
    colvals[3] = 4;
    colvals[4] = 5;
    colvals[5] = 6;
    colvals[6] = 8;
    colvals[7] = 12;
    colvals[8] = 0;
    colvals[9] = 1;
    colvals[10] = 3;
    colvals[11] = 4;
    colvals[12] = 6;
    colvals[13] = 8;
    colvals[14] = 12;
    colvals[15] = 3;
    colvals[16] = 4;
    colvals[17] = 10;
    colvals[18] = 11;
    colvals[19] = 12;
    colvals[20] = 0;
    colvals[21] = 1;
    colvals[22] = 3;
    colvals[23] = 4;
    colvals[24] = 5;
    colvals[25] = 6;
    colvals[26] = 7;
    colvals[27] = 8;
    colvals[28] = 12;
    colvals[29] = 0;
    colvals[30] = 1;
    colvals[31] = 3;
    colvals[32] = 4;
    colvals[33] = 5;
    colvals[34] = 6;
    colvals[35] = 7;
    colvals[36] = 8;
    colvals[37] = 12;
    colvals[38] = 0;
    colvals[39] = 3;
    colvals[40] = 4;
    colvals[41] = 5;
    colvals[42] = 7;
    colvals[43] = 12;
    colvals[44] = 0;
    colvals[45] = 1;
    colvals[46] = 3;
    colvals[47] = 4;
    colvals[48] = 5;
    colvals[49] = 6;
    colvals[50] = 7;
    colvals[51] = 8;
    colvals[52] = 12;
    colvals[53] = 3;
    colvals[54] = 4;
    colvals[55] = 5;
    colvals[56] = 6;
    colvals[57] = 7;
    colvals[58] = 12;
    colvals[59] = 0;
    colvals[60] = 1;
    colvals[61] = 3;
    colvals[62] = 4;
    colvals[63] = 5;
    colvals[64] = 6;
    colvals[65] = 8;
    colvals[66] = 9;
    colvals[67] = 10;
    colvals[68] = 12;
    colvals[69] = 9;
    colvals[70] = 10;
    colvals[71] = 11;
    colvals[72] = 12;
    colvals[73] = 10;
    colvals[74] = 11;
    colvals[75] = 12;
    colvals[76] = 0;
    colvals[77] = 1;
    colvals[78] = 3;
    colvals[79] = 4;
    colvals[80] = 5;
    colvals[81] = 6;
    colvals[82] = 7;
    colvals[83] = 9;
    colvals[84] = 10;
    colvals[85] = 11;
    colvals[86] = 12;
    colvals[87] = 3;
    colvals[88] = 4;
    colvals[89] = 9;
    colvals[90] = 10;
    colvals[91] = 11;
    colvals[92] = 12;
    
    // value of each non-zero element
    data[0] = 0.0 - k[29]*y[IDX_HII] - k[33]*y[IDX_H2I] - k[34]*y[IDX_H2I] -
        k[36]*y[IDX_HM];
    data[1] = 0.0 + k[30]*y[IDX_HI] + k[37]*y[IDX_eM];
    data[2] = 0.0 + k[30]*y[IDX_DII] + k[35]*y[IDX_HDI];
    data[3] = 0.0 - k[29]*y[IDX_DI];
    data[4] = 0.0 - k[36]*y[IDX_DI];
    data[5] = 0.0 - k[33]*y[IDX_DI] - k[34]*y[IDX_DI];
    data[6] = 0.0 + k[35]*y[IDX_HI];
    data[7] = 0.0 + k[37]*y[IDX_DII];
    data[8] = 0.0 + k[29]*y[IDX_HII];
    data[9] = 0.0 - k[30]*y[IDX_HI] - k[31]*y[IDX_H2I] - k[37]*y[IDX_eM];
    data[10] = 0.0 - k[30]*y[IDX_DII];
    data[11] = 0.0 + k[29]*y[IDX_DI] + k[32]*y[IDX_HDI];
    data[12] = 0.0 - k[31]*y[IDX_DII];
    data[13] = 0.0 + k[32]*y[IDX_HII];
    data[14] = 0.0 - k[37]*y[IDX_DII];
    data[15] = 0.0 + k[8]*y[IDX_eM] + k[11]*y[IDX_HII] + k[12]*y[IDX_HII];
    data[16] = 0.0 + k[1]*y[IDX_eM] + k[2]*y[IDX_eM] + k[11]*y[IDX_HI] +
        k[12]*y[IDX_HI];
    data[17] = 0.0 + k[4]*y[IDX_eM] + k[5]*y[IDX_eM];
    data[18] = 0.0 + k[7]*y[IDX_eM];
    data[19] = 0.0 + k[1]*y[IDX_HII] + k[2]*y[IDX_HII] + k[4]*y[IDX_HeII] +
        k[5]*y[IDX_HeII] + k[7]*y[IDX_HeIII] + k[8]*y[IDX_HI];
    data[20] = 0.0 + k[29]*y[IDX_HII] + k[33]*y[IDX_H2I] + k[34]*y[IDX_H2I];
    data[21] = 0.0 - k[30]*y[IDX_HI];
    data[22] = 0.0 - k[0]*y[IDX_eM] - k[8]*y[IDX_eM] - k[9]*y[IDX_HM] -
        k[10]*y[IDX_HM] - k[11]*y[IDX_HII] - k[12]*y[IDX_HII] -
        k[13]*y[IDX_H2II] - k[16]*y[IDX_H2I] + k[16]*y[IDX_H2I] +
        k[16]*y[IDX_H2I] + k[16]*y[IDX_H2I] - k[18]*y[IDX_HM] + k[18]*y[IDX_HM]
        + k[18]*y[IDX_HM] - k[19]*y[IDX_HM] + k[19]*y[IDX_HM] + k[19]*y[IDX_HM]
        - k[25]*y[IDX_HI]*y[IDX_HI] - k[25]*y[IDX_HI]*y[IDX_HI] -
        k[25]*y[IDX_HI]*y[IDX_HI] - k[25]*y[IDX_HI]*y[IDX_HI] -
        k[25]*y[IDX_HI]*y[IDX_HI] - k[25]*y[IDX_HI]*y[IDX_HI] -
        k[25]*y[IDX_HI]*y[IDX_HI] - k[25]*y[IDX_HI]*y[IDX_HI] -
        k[25]*y[IDX_HI]*y[IDX_HI] + k[25]*y[IDX_HI]*y[IDX_HI] +
        k[25]*y[IDX_HI]*y[IDX_HI] + k[25]*y[IDX_HI]*y[IDX_HI] -
        k[26]*y[IDX_HI]*y[IDX_HI] - k[26]*y[IDX_HI]*y[IDX_HI] -
        k[26]*y[IDX_HI]*y[IDX_HI] - k[26]*y[IDX_HI]*y[IDX_HI] -
        k[26]*y[IDX_HI]*y[IDX_HI] - k[26]*y[IDX_HI]*y[IDX_HI] -
        k[26]*y[IDX_HI]*y[IDX_HI] - k[26]*y[IDX_HI]*y[IDX_HI] -
        k[26]*y[IDX_HI]*y[IDX_HI] + k[26]*y[IDX_HI]*y[IDX_HI] +
        k[26]*y[IDX_HI]*y[IDX_HI] + k[26]*y[IDX_HI]*y[IDX_HI] -
        k[27]*y[IDX_H2I]*y[IDX_HI] - k[27]*y[IDX_H2I]*y[IDX_HI] -
        k[27]*y[IDX_H2I]*y[IDX_HI] - k[27]*y[IDX_H2I]*y[IDX_HI] -
        k[28]*y[IDX_H2I]*y[IDX_HI] - k[28]*y[IDX_H2I]*y[IDX_HI] -
        k[28]*y[IDX_H2I]*y[IDX_HI] - k[28]*y[IDX_H2I]*y[IDX_HI] -
        k[30]*y[IDX_DII] - k[35]*y[IDX_HDI];
    data[23] = 0.0 + k[1]*y[IDX_eM] + k[2]*y[IDX_eM] - k[11]*y[IDX_HI] -
        k[12]*y[IDX_HI] + k[14]*y[IDX_H2I] + k[20]*y[IDX_HM] + k[20]*y[IDX_HM] +
        k[29]*y[IDX_DI];
    data[24] = 0.0 - k[9]*y[IDX_HI] - k[10]*y[IDX_HI] + k[17]*y[IDX_eM] -
        k[18]*y[IDX_HI] + k[18]*y[IDX_HI] + k[18]*y[IDX_HI] - k[19]*y[IDX_HI] +
        k[19]*y[IDX_HI] + k[19]*y[IDX_HI] + k[20]*y[IDX_HII] + k[20]*y[IDX_HII]
        + k[24]*y[IDX_H2II];
    data[25] = 0.0 + k[14]*y[IDX_HII] + k[15]*y[IDX_eM] + k[15]*y[IDX_eM] -
        k[16]*y[IDX_HI] + k[16]*y[IDX_HI] + k[16]*y[IDX_HI] + k[16]*y[IDX_HI] -
        k[27]*y[IDX_HI]*y[IDX_HI] - k[27]*y[IDX_HI]*y[IDX_HI] -
        k[28]*y[IDX_HI]*y[IDX_HI] - k[28]*y[IDX_HI]*y[IDX_HI] + k[33]*y[IDX_DI]
        + k[34]*y[IDX_DI];
    data[26] = 0.0 - k[13]*y[IDX_HI] + k[22]*y[IDX_eM] + k[22]*y[IDX_eM] +
        k[23]*y[IDX_eM] + k[23]*y[IDX_eM] + k[24]*y[IDX_HM];
    data[27] = 0.0 - k[35]*y[IDX_HI];
    data[28] = 0.0 - k[0]*y[IDX_HI] + k[1]*y[IDX_HII] + k[2]*y[IDX_HII] -
        k[8]*y[IDX_HI] + k[15]*y[IDX_H2I] + k[15]*y[IDX_H2I] + k[17]*y[IDX_HM] +
        k[22]*y[IDX_H2II] + k[22]*y[IDX_H2II] + k[23]*y[IDX_H2II] +
        k[23]*y[IDX_H2II];
    data[29] = 0.0 - k[29]*y[IDX_HII];
    data[30] = 0.0 + k[30]*y[IDX_HI] + k[31]*y[IDX_H2I];
    data[31] = 0.0 + k[0]*y[IDX_eM] - k[11]*y[IDX_HII] - k[12]*y[IDX_HII] +
        k[13]*y[IDX_H2II] + k[30]*y[IDX_DII];
    data[32] = 0.0 - k[1]*y[IDX_eM] - k[2]*y[IDX_eM] - k[11]*y[IDX_HI] -
        k[12]*y[IDX_HI] - k[14]*y[IDX_H2I] - k[20]*y[IDX_HM] - k[21]*y[IDX_HM] -
        k[29]*y[IDX_DI] - k[32]*y[IDX_HDI];
    data[33] = 0.0 - k[20]*y[IDX_HII] - k[21]*y[IDX_HII];
    data[34] = 0.0 - k[14]*y[IDX_HII] + k[31]*y[IDX_DII];
    data[35] = 0.0 + k[13]*y[IDX_HI];
    data[36] = 0.0 - k[32]*y[IDX_HII];
    data[37] = 0.0 + k[0]*y[IDX_HI] - k[1]*y[IDX_HII] - k[2]*y[IDX_HII];
    data[38] = 0.0 - k[36]*y[IDX_HM];
    data[39] = 0.0 + k[8]*y[IDX_eM] - k[9]*y[IDX_HM] - k[10]*y[IDX_HM] -
        k[18]*y[IDX_HM] - k[19]*y[IDX_HM];
    data[40] = 0.0 - k[20]*y[IDX_HM] - k[21]*y[IDX_HM];
    data[41] = 0.0 - k[9]*y[IDX_HI] - k[10]*y[IDX_HI] - k[17]*y[IDX_eM] -
        k[18]*y[IDX_HI] - k[19]*y[IDX_HI] - k[20]*y[IDX_HII] - k[21]*y[IDX_HII]
        - k[24]*y[IDX_H2II] - k[36]*y[IDX_DI];
    data[42] = 0.0 - k[24]*y[IDX_HM];
    data[43] = 0.0 + k[8]*y[IDX_HI] - k[17]*y[IDX_HM];
    data[44] = 0.0 - k[33]*y[IDX_H2I] - k[34]*y[IDX_H2I];
    data[45] = 0.0 - k[31]*y[IDX_H2I];
    data[46] = 0.0 + k[9]*y[IDX_HM] + k[10]*y[IDX_HM] + k[13]*y[IDX_H2II] -
        k[16]*y[IDX_H2I] + k[25]*y[IDX_HI]*y[IDX_HI] + k[25]*y[IDX_HI]*y[IDX_HI]
        + k[25]*y[IDX_HI]*y[IDX_HI] + k[26]*y[IDX_HI]*y[IDX_HI] +
        k[26]*y[IDX_HI]*y[IDX_HI] + k[26]*y[IDX_HI]*y[IDX_HI] -
        k[27]*y[IDX_H2I]*y[IDX_HI] - k[27]*y[IDX_H2I]*y[IDX_HI] +
        k[27]*y[IDX_H2I]*y[IDX_HI] + k[27]*y[IDX_H2I]*y[IDX_HI] +
        k[27]*y[IDX_H2I]*y[IDX_HI] + k[27]*y[IDX_H2I]*y[IDX_HI] -
        k[28]*y[IDX_H2I]*y[IDX_HI] - k[28]*y[IDX_H2I]*y[IDX_HI] +
        k[28]*y[IDX_H2I]*y[IDX_HI] + k[28]*y[IDX_H2I]*y[IDX_HI] +
        k[28]*y[IDX_H2I]*y[IDX_HI] + k[28]*y[IDX_H2I]*y[IDX_HI] +
        k[35]*y[IDX_HDI];
    data[47] = 0.0 - k[14]*y[IDX_H2I] + k[32]*y[IDX_HDI];
    data[48] = 0.0 + k[9]*y[IDX_HI] + k[10]*y[IDX_HI] + k[24]*y[IDX_H2II];
    data[49] = 0.0 - k[14]*y[IDX_HII] - k[15]*y[IDX_eM] - k[16]*y[IDX_HI] -
        k[27]*y[IDX_HI]*y[IDX_HI] + k[27]*y[IDX_HI]*y[IDX_HI] +
        k[27]*y[IDX_HI]*y[IDX_HI] - k[28]*y[IDX_HI]*y[IDX_HI] +
        k[28]*y[IDX_HI]*y[IDX_HI] + k[28]*y[IDX_HI]*y[IDX_HI] - k[31]*y[IDX_DII]
        - k[33]*y[IDX_DI] - k[34]*y[IDX_DI];
    data[50] = 0.0 + k[13]*y[IDX_HI] + k[24]*y[IDX_HM];
    data[51] = 0.0 + k[32]*y[IDX_HII] + k[35]*y[IDX_HI];
    data[52] = 0.0 - k[15]*y[IDX_H2I];
    data[53] = 0.0 + k[11]*y[IDX_HII] + k[12]*y[IDX_HII] -
        k[13]*y[IDX_H2II];
    data[54] = 0.0 + k[11]*y[IDX_HI] + k[12]*y[IDX_HI] + k[14]*y[IDX_H2I] +
        k[21]*y[IDX_HM];
    data[55] = 0.0 + k[21]*y[IDX_HII] - k[24]*y[IDX_H2II];
    data[56] = 0.0 + k[14]*y[IDX_HII];
    data[57] = 0.0 - k[13]*y[IDX_HI] - k[22]*y[IDX_eM] - k[23]*y[IDX_eM] -
        k[24]*y[IDX_HM];
    data[58] = 0.0 - k[22]*y[IDX_H2II] - k[23]*y[IDX_H2II];
    data[59] = 0.0 + k[33]*y[IDX_H2I] + k[34]*y[IDX_H2I] + k[36]*y[IDX_HM];
    data[60] = 0.0 + k[31]*y[IDX_H2I];
    data[61] = 0.0 - k[35]*y[IDX_HDI];
    data[62] = 0.0 - k[32]*y[IDX_HDI];
    data[63] = 0.0 + k[36]*y[IDX_DI];
    data[64] = 0.0 + k[31]*y[IDX_DII] + k[33]*y[IDX_DI] + k[34]*y[IDX_DI];
    data[65] = 0.0 - k[32]*y[IDX_HII] - k[35]*y[IDX_HI];
    data[66] = 0.0 - k[3]*y[IDX_eM];
    data[67] = 0.0 + k[4]*y[IDX_eM] + k[5]*y[IDX_eM];
    data[68] = 0.0 - k[3]*y[IDX_HeI] + k[4]*y[IDX_HeII] + k[5]*y[IDX_HeII];
    data[69] = 0.0 + k[3]*y[IDX_eM];
    data[70] = 0.0 - k[4]*y[IDX_eM] - k[5]*y[IDX_eM] - k[6]*y[IDX_eM];
    data[71] = 0.0 + k[7]*y[IDX_eM];
    data[72] = 0.0 + k[3]*y[IDX_HeI] - k[4]*y[IDX_HeII] - k[5]*y[IDX_HeII] -
        k[6]*y[IDX_HeII] + k[7]*y[IDX_HeIII];
    data[73] = 0.0 + k[6]*y[IDX_eM];
    data[74] = 0.0 - k[7]*y[IDX_eM];
    data[75] = 0.0 + k[6]*y[IDX_HeII] - k[7]*y[IDX_HeIII];
    data[76] = 0.0 + k[36]*y[IDX_HM];
    data[77] = 0.0 - k[37]*y[IDX_eM];
    data[78] = 0.0 - k[0]*y[IDX_eM] + k[0]*y[IDX_eM] + k[0]*y[IDX_eM] -
        k[8]*y[IDX_eM] + k[9]*y[IDX_HM] + k[10]*y[IDX_HM] + k[18]*y[IDX_HM] +
        k[19]*y[IDX_HM];
    data[79] = 0.0 - k[1]*y[IDX_eM] - k[2]*y[IDX_eM] + k[21]*y[IDX_HM];
    data[80] = 0.0 + k[9]*y[IDX_HI] + k[10]*y[IDX_HI] - k[17]*y[IDX_eM] +
        k[17]*y[IDX_eM] + k[17]*y[IDX_eM] + k[18]*y[IDX_HI] + k[19]*y[IDX_HI] +
        k[21]*y[IDX_HII] + k[36]*y[IDX_DI];
    data[81] = 0.0 - k[15]*y[IDX_eM] + k[15]*y[IDX_eM];
    data[82] = 0.0 - k[22]*y[IDX_eM] - k[23]*y[IDX_eM];
    data[83] = 0.0 - k[3]*y[IDX_eM] + k[3]*y[IDX_eM] + k[3]*y[IDX_eM];
    data[84] = 0.0 - k[4]*y[IDX_eM] - k[5]*y[IDX_eM] - k[6]*y[IDX_eM] +
        k[6]*y[IDX_eM] + k[6]*y[IDX_eM];
    data[85] = 0.0 - k[7]*y[IDX_eM];
    data[86] = 0.0 - k[0]*y[IDX_HI] + k[0]*y[IDX_HI] + k[0]*y[IDX_HI] -
        k[1]*y[IDX_HII] - k[2]*y[IDX_HII] - k[3]*y[IDX_HeI] + k[3]*y[IDX_HeI] +
        k[3]*y[IDX_HeI] - k[4]*y[IDX_HeII] - k[5]*y[IDX_HeII] - k[6]*y[IDX_HeII]
        + k[6]*y[IDX_HeII] + k[6]*y[IDX_HeII] - k[7]*y[IDX_HeIII] -
        k[8]*y[IDX_HI] - k[15]*y[IDX_H2I] + k[15]*y[IDX_H2I] - k[17]*y[IDX_HM] +
        k[17]*y[IDX_HM] + k[17]*y[IDX_HM] - k[22]*y[IDX_H2II] -
        k[23]*y[IDX_H2II] - k[37]*y[IDX_DII];
    data[87] = 0.0 - kc[0]*y[IDX_eM] - kc[8]*y[IDX_eM]*y[IDX_eM];
    data[88] = 0.0 - kc[4]*y[IDX_eM];
    data[89] = 0.0 - kc[1]*y[IDX_eM];
    data[90] = 0.0 - kc[2]*y[IDX_eM] - kc[3]*y[IDX_eM]*y[IDX_eM] -
        kc[5]*y[IDX_eM] - kc[6]*y[IDX_eM] - kc[9]*y[IDX_eM] - kc[10]*y[IDX_eM];
    data[91] = 0.0 - kc[7]*y[IDX_eM];
    data[92] = (gamma - 1.0) * (0.0 - kc[0]*y[IDX_HI] - kc[1]*y[IDX_HeI] -
        kc[2]*y[IDX_HeII] - kc[3]*y[IDX_HeII]*y[IDX_eM] -
        kc[3]*y[IDX_HeII]*y[IDX_eM] - kc[4]*y[IDX_HII] - kc[5]*y[IDX_HeII] -
        kc[6]*y[IDX_HeII] - kc[7]*y[IDX_HeIII] - kc[8]*y[IDX_HI]*y[IDX_eM] -
        kc[8]*y[IDX_HI]*y[IDX_eM] - kc[9]*y[IDX_HeII] - kc[10]*y[IDX_HeII] ) /
        kerg / GetNumDens(y);
    
    // clang-format on

    /* */

    return NAUNET_SUCCESS;
}