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

int Jac(realtype t, N_Vector u, N_Vector fu, SUNMatrix jmatrix, void *user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    /* */
    realtype *y            = N_VGetArrayPointer(u);
    sunindextype *rowptrs  = SUNSparseMatrix_IndexPointers(jmatrix);
    sunindextype *colvals  = SUNSparseMatrix_IndexValues(jmatrix);
    realtype *data         = SUNSparseMatrix_Data(jmatrix);
    NaunetData *u_data     = (NaunetData *)user_data;

    // clang-format off
    realtype nH = u_data->nH;
    realtype Tgas = u_data->Tgas;
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
    // number of non-zero elements in each row
    rowptrs[0] = 0;
    rowptrs[1] = 8;
    rowptrs[2] = 15;
    rowptrs[3] = 24;
    rowptrs[4] = 33;
    rowptrs[5] = 39;
    rowptrs[6] = 48;
    rowptrs[7] = 54;
    rowptrs[8] = 61;
    rowptrs[9] = 64;
    rowptrs[10] = 68;
    rowptrs[11] = 71;
    rowptrs[12] = 82;
    rowptrs[13] = 88;
    
    // the column index of non-zero elements
    colvals[0] = 0;
    colvals[1] = 1;
    colvals[2] = 2;
    colvals[3] = 3;
    colvals[4] = 4;
    colvals[5] = 5;
    colvals[6] = 7;
    colvals[7] = 11;
    colvals[8] = 0;
    colvals[9] = 1;
    colvals[10] = 2;
    colvals[11] = 3;
    colvals[12] = 5;
    colvals[13] = 7;
    colvals[14] = 11;
    colvals[15] = 0;
    colvals[16] = 1;
    colvals[17] = 2;
    colvals[18] = 3;
    colvals[19] = 4;
    colvals[20] = 5;
    colvals[21] = 6;
    colvals[22] = 7;
    colvals[23] = 11;
    colvals[24] = 0;
    colvals[25] = 1;
    colvals[26] = 2;
    colvals[27] = 3;
    colvals[28] = 4;
    colvals[29] = 5;
    colvals[30] = 6;
    colvals[31] = 7;
    colvals[32] = 11;
    colvals[33] = 0;
    colvals[34] = 2;
    colvals[35] = 3;
    colvals[36] = 4;
    colvals[37] = 6;
    colvals[38] = 11;
    colvals[39] = 0;
    colvals[40] = 1;
    colvals[41] = 2;
    colvals[42] = 3;
    colvals[43] = 4;
    colvals[44] = 5;
    colvals[45] = 6;
    colvals[46] = 7;
    colvals[47] = 11;
    colvals[48] = 2;
    colvals[49] = 3;
    colvals[50] = 4;
    colvals[51] = 5;
    colvals[52] = 6;
    colvals[53] = 11;
    colvals[54] = 0;
    colvals[55] = 1;
    colvals[56] = 2;
    colvals[57] = 3;
    colvals[58] = 4;
    colvals[59] = 5;
    colvals[60] = 7;
    colvals[61] = 8;
    colvals[62] = 9;
    colvals[63] = 11;
    colvals[64] = 8;
    colvals[65] = 9;
    colvals[66] = 10;
    colvals[67] = 11;
    colvals[68] = 9;
    colvals[69] = 10;
    colvals[70] = 11;
    colvals[71] = 0;
    colvals[72] = 1;
    colvals[73] = 2;
    colvals[74] = 3;
    colvals[75] = 4;
    colvals[76] = 5;
    colvals[77] = 6;
    colvals[78] = 8;
    colvals[79] = 9;
    colvals[80] = 10;
    colvals[81] = 11;
    colvals[82] = 2;
    colvals[83] = 3;
    colvals[84] = 8;
    colvals[85] = 9;
    colvals[86] = 10;
    colvals[87] = 11;
    
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
    data[15] = 0.0 + k[29]*y[IDX_HII] + k[33]*y[IDX_H2I] + k[34]*y[IDX_H2I];
    data[16] = 0.0 - k[30]*y[IDX_HI];
    data[17] = 0.0 - k[0]*y[IDX_eM] - k[8]*y[IDX_eM] - k[9]*y[IDX_HM] - k[10]*y[IDX_HM]
        - k[11]*y[IDX_HII] - k[12]*y[IDX_HII] - k[13]*y[IDX_H2II] -
        k[16]*y[IDX_H2I] + k[16]*y[IDX_H2I] + k[16]*y[IDX_H2I] +
        k[16]*y[IDX_H2I] - k[18]*y[IDX_HM] + k[18]*y[IDX_HM] + k[18]*y[IDX_HM] -
        k[19]*y[IDX_HM] + k[19]*y[IDX_HM] + k[19]*y[IDX_HM] -
        k[25]*y[IDX_HI]*y[IDX_HI] - k[25]*y[IDX_HI]*y[IDX_HI] -
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
    data[18] = 0.0 + k[1]*y[IDX_eM] + k[2]*y[IDX_eM] - k[11]*y[IDX_HI] -
        k[12]*y[IDX_HI] + k[14]*y[IDX_H2I] + k[20]*y[IDX_HM] + k[20]*y[IDX_HM] +
        k[29]*y[IDX_DI];
    data[19] = 0.0 - k[9]*y[IDX_HI] - k[10]*y[IDX_HI] + k[17]*y[IDX_eM] -
        k[18]*y[IDX_HI] + k[18]*y[IDX_HI] + k[18]*y[IDX_HI] - k[19]*y[IDX_HI] +
        k[19]*y[IDX_HI] + k[19]*y[IDX_HI] + k[20]*y[IDX_HII] + k[20]*y[IDX_HII]
        + k[24]*y[IDX_H2II];
    data[20] = 0.0 + k[14]*y[IDX_HII] + k[15]*y[IDX_eM] + k[15]*y[IDX_eM] -
        k[16]*y[IDX_HI] + k[16]*y[IDX_HI] + k[16]*y[IDX_HI] + k[16]*y[IDX_HI] -
        k[27]*y[IDX_HI]*y[IDX_HI] - k[27]*y[IDX_HI]*y[IDX_HI] -
        k[28]*y[IDX_HI]*y[IDX_HI] - k[28]*y[IDX_HI]*y[IDX_HI] + k[33]*y[IDX_DI]
        + k[34]*y[IDX_DI];
    data[21] = 0.0 - k[13]*y[IDX_HI] + k[22]*y[IDX_eM] + k[22]*y[IDX_eM] +
        k[23]*y[IDX_eM] + k[23]*y[IDX_eM] + k[24]*y[IDX_HM];
    data[22] = 0.0 - k[35]*y[IDX_HI];
    data[23] = 0.0 - k[0]*y[IDX_HI] + k[1]*y[IDX_HII] + k[2]*y[IDX_HII] -
        k[8]*y[IDX_HI] + k[15]*y[IDX_H2I] + k[15]*y[IDX_H2I] + k[17]*y[IDX_HM] +
        k[22]*y[IDX_H2II] + k[22]*y[IDX_H2II] + k[23]*y[IDX_H2II] +
        k[23]*y[IDX_H2II];
    data[24] = 0.0 - k[29]*y[IDX_HII];
    data[25] = 0.0 + k[30]*y[IDX_HI] + k[31]*y[IDX_H2I];
    data[26] = 0.0 + k[0]*y[IDX_eM] - k[11]*y[IDX_HII] - k[12]*y[IDX_HII] +
        k[13]*y[IDX_H2II] + k[30]*y[IDX_DII];
    data[27] = 0.0 - k[1]*y[IDX_eM] - k[2]*y[IDX_eM] - k[11]*y[IDX_HI] -
        k[12]*y[IDX_HI] - k[14]*y[IDX_H2I] - k[20]*y[IDX_HM] - k[21]*y[IDX_HM] -
        k[29]*y[IDX_DI] - k[32]*y[IDX_HDI];
    data[28] = 0.0 - k[20]*y[IDX_HII] - k[21]*y[IDX_HII];
    data[29] = 0.0 - k[14]*y[IDX_HII] + k[31]*y[IDX_DII];
    data[30] = 0.0 + k[13]*y[IDX_HI];
    data[31] = 0.0 - k[32]*y[IDX_HII];
    data[32] = 0.0 + k[0]*y[IDX_HI] - k[1]*y[IDX_HII] - k[2]*y[IDX_HII];
    data[33] = 0.0 - k[36]*y[IDX_HM];
    data[34] = 0.0 + k[8]*y[IDX_eM] - k[9]*y[IDX_HM] - k[10]*y[IDX_HM] -
        k[18]*y[IDX_HM] - k[19]*y[IDX_HM];
    data[35] = 0.0 - k[20]*y[IDX_HM] - k[21]*y[IDX_HM];
    data[36] = 0.0 - k[9]*y[IDX_HI] - k[10]*y[IDX_HI] - k[17]*y[IDX_eM] -
        k[18]*y[IDX_HI] - k[19]*y[IDX_HI] - k[20]*y[IDX_HII] - k[21]*y[IDX_HII]
        - k[24]*y[IDX_H2II] - k[36]*y[IDX_DI];
    data[37] = 0.0 - k[24]*y[IDX_HM];
    data[38] = 0.0 + k[8]*y[IDX_HI] - k[17]*y[IDX_HM];
    data[39] = 0.0 - k[33]*y[IDX_H2I] - k[34]*y[IDX_H2I];
    data[40] = 0.0 - k[31]*y[IDX_H2I];
    data[41] = 0.0 + k[9]*y[IDX_HM] + k[10]*y[IDX_HM] + k[13]*y[IDX_H2II] -
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
    data[42] = 0.0 - k[14]*y[IDX_H2I] + k[32]*y[IDX_HDI];
    data[43] = 0.0 + k[9]*y[IDX_HI] + k[10]*y[IDX_HI] + k[24]*y[IDX_H2II];
    data[44] = 0.0 - k[14]*y[IDX_HII] - k[15]*y[IDX_eM] - k[16]*y[IDX_HI] -
        k[27]*y[IDX_HI]*y[IDX_HI] + k[27]*y[IDX_HI]*y[IDX_HI] +
        k[27]*y[IDX_HI]*y[IDX_HI] - k[28]*y[IDX_HI]*y[IDX_HI] +
        k[28]*y[IDX_HI]*y[IDX_HI] + k[28]*y[IDX_HI]*y[IDX_HI] - k[31]*y[IDX_DII]
        - k[33]*y[IDX_DI] - k[34]*y[IDX_DI];
    data[45] = 0.0 + k[13]*y[IDX_HI] + k[24]*y[IDX_HM];
    data[46] = 0.0 + k[32]*y[IDX_HII] + k[35]*y[IDX_HI];
    data[47] = 0.0 - k[15]*y[IDX_H2I];
    data[48] = 0.0 + k[11]*y[IDX_HII] + k[12]*y[IDX_HII] - k[13]*y[IDX_H2II];
    data[49] = 0.0 + k[11]*y[IDX_HI] + k[12]*y[IDX_HI] + k[14]*y[IDX_H2I] +
        k[21]*y[IDX_HM];
    data[50] = 0.0 + k[21]*y[IDX_HII] - k[24]*y[IDX_H2II];
    data[51] = 0.0 + k[14]*y[IDX_HII];
    data[52] = 0.0 - k[13]*y[IDX_HI] - k[22]*y[IDX_eM] - k[23]*y[IDX_eM] -
        k[24]*y[IDX_HM];
    data[53] = 0.0 - k[22]*y[IDX_H2II] - k[23]*y[IDX_H2II];
    data[54] = 0.0 + k[33]*y[IDX_H2I] + k[34]*y[IDX_H2I] + k[36]*y[IDX_HM];
    data[55] = 0.0 + k[31]*y[IDX_H2I];
    data[56] = 0.0 - k[35]*y[IDX_HDI];
    data[57] = 0.0 - k[32]*y[IDX_HDI];
    data[58] = 0.0 + k[36]*y[IDX_DI];
    data[59] = 0.0 + k[31]*y[IDX_DII] + k[33]*y[IDX_DI] + k[34]*y[IDX_DI];
    data[60] = 0.0 - k[32]*y[IDX_HII] - k[35]*y[IDX_HI];
    data[61] = 0.0 - k[3]*y[IDX_eM];
    data[62] = 0.0 + k[4]*y[IDX_eM] + k[5]*y[IDX_eM];
    data[63] = 0.0 - k[3]*y[IDX_HeI] + k[4]*y[IDX_HeII] + k[5]*y[IDX_HeII];
    data[64] = 0.0 + k[3]*y[IDX_eM];
    data[65] = 0.0 - k[4]*y[IDX_eM] - k[5]*y[IDX_eM] - k[6]*y[IDX_eM];
    data[66] = 0.0 + k[7]*y[IDX_eM];
    data[67] = 0.0 + k[3]*y[IDX_HeI] - k[4]*y[IDX_HeII] - k[5]*y[IDX_HeII] -
        k[6]*y[IDX_HeII] + k[7]*y[IDX_HeIII];
    data[68] = 0.0 + k[6]*y[IDX_eM];
    data[69] = 0.0 - k[7]*y[IDX_eM];
    data[70] = 0.0 + k[6]*y[IDX_HeII] - k[7]*y[IDX_HeIII];
    data[71] = 0.0 + k[36]*y[IDX_HM];
    data[72] = 0.0 - k[37]*y[IDX_eM];
    data[73] = 0.0 - k[0]*y[IDX_eM] + k[0]*y[IDX_eM] + k[0]*y[IDX_eM] - k[8]*y[IDX_eM]
        + k[9]*y[IDX_HM] + k[10]*y[IDX_HM] + k[18]*y[IDX_HM] + k[19]*y[IDX_HM];
    data[74] = 0.0 - k[1]*y[IDX_eM] - k[2]*y[IDX_eM] + k[21]*y[IDX_HM];
    data[75] = 0.0 + k[9]*y[IDX_HI] + k[10]*y[IDX_HI] - k[17]*y[IDX_eM] +
        k[17]*y[IDX_eM] + k[17]*y[IDX_eM] + k[18]*y[IDX_HI] + k[19]*y[IDX_HI] +
        k[21]*y[IDX_HII] + k[36]*y[IDX_DI];
    data[76] = 0.0 - k[15]*y[IDX_eM] + k[15]*y[IDX_eM];
    data[77] = 0.0 - k[22]*y[IDX_eM] - k[23]*y[IDX_eM];
    data[78] = 0.0 - k[3]*y[IDX_eM] + k[3]*y[IDX_eM] + k[3]*y[IDX_eM];
    data[79] = 0.0 - k[4]*y[IDX_eM] - k[5]*y[IDX_eM] - k[6]*y[IDX_eM] + k[6]*y[IDX_eM]
        + k[6]*y[IDX_eM];
    data[80] = 0.0 - k[7]*y[IDX_eM];
    data[81] = 0.0 - k[0]*y[IDX_HI] + k[0]*y[IDX_HI] + k[0]*y[IDX_HI] - k[1]*y[IDX_HII]
        - k[2]*y[IDX_HII] - k[3]*y[IDX_HeI] + k[3]*y[IDX_HeI] + k[3]*y[IDX_HeI]
        - k[4]*y[IDX_HeII] - k[5]*y[IDX_HeII] - k[6]*y[IDX_HeII] +
        k[6]*y[IDX_HeII] + k[6]*y[IDX_HeII] - k[7]*y[IDX_HeIII] - k[8]*y[IDX_HI]
        - k[15]*y[IDX_H2I] + k[15]*y[IDX_H2I] - k[17]*y[IDX_HM] +
        k[17]*y[IDX_HM] + k[17]*y[IDX_HM] - k[22]*y[IDX_H2II] -
        k[23]*y[IDX_H2II] - k[37]*y[IDX_DII];
    data[82] = (gamma - 1.0) * (0.0 - kc[0]*y[IDX_eM] - kc[8]*y[IDX_eM]*y[IDX_eM] ) /
        kerg / GetNumDens(y);
    data[83] = (gamma - 1.0) * (0.0 - kc[4]*y[IDX_eM] ) / kerg / GetNumDens(y);
    data[84] = (gamma - 1.0) * (0.0 - kc[1]*y[IDX_eM] ) / kerg / GetNumDens(y);
    data[85] = (gamma - 1.0) * (0.0 - kc[2]*y[IDX_eM] - kc[3]*y[IDX_eM]*y[IDX_eM] -
        kc[5]*y[IDX_eM] - kc[6]*y[IDX_eM] - kc[9]*y[IDX_eM] - kc[10]*y[IDX_eM] )
        / kerg / GetNumDens(y);
    data[86] = (gamma - 1.0) * (0.0 - kc[7]*y[IDX_eM] ) / kerg / GetNumDens(y);
    data[87] = (gamma - 1.0) * (0.0 - kc[0]*y[IDX_HI] - kc[1]*y[IDX_HeI] -
        kc[2]*y[IDX_HeII] - kc[3]*y[IDX_HeII]*y[IDX_eM] -
        kc[3]*y[IDX_HeII]*y[IDX_eM] - kc[4]*y[IDX_HII] - kc[5]*y[IDX_HeII] -
        kc[6]*y[IDX_HeII] - kc[7]*y[IDX_HeIII] - kc[8]*y[IDX_HI]*y[IDX_eM] -
        kc[8]*y[IDX_HI]*y[IDX_eM] - kc[9]*y[IDX_HeII] - kc[10]*y[IDX_HeII] ) /
        kerg / GetNumDens(y);
    
    // clang-format on

    /* */

    return NAUNET_SUCCESS;
}