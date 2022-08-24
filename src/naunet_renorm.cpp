#include "naunet_macros.h"
#include "naunet_physics.h"
#include "naunet_renorm.h"

// clang-format off
int InitRenorm(realtype *ab, SUNMatrix A) {
    // clang-format on
    realtype Hnuclei = GetHNuclei(ab);

    // clang-format off
            
    IJth(A, IDX_ELEM_He, IDX_ELEM_He) = 0.0 + 4.0 * ab[IDX_HeI] / 4.0 / Hnuclei +
                                    4.0 * ab[IDX_HeIII] / 4.0 / Hnuclei + 4.0 *
                                    ab[IDX_HeII] / 4.0 / Hnuclei;
    IJth(A, IDX_ELEM_He, IDX_ELEM_D) = 0.0;
    IJth(A, IDX_ELEM_He, IDX_ELEM_H) = 0.0;
    IJth(A, IDX_ELEM_D, IDX_ELEM_He) = 0.0;
    IJth(A, IDX_ELEM_D, IDX_ELEM_D) = 0.0 + 2.0 * ab[IDX_DII] / 2.0 / Hnuclei +
                                    2.0 * ab[IDX_DI] / 2.0 / Hnuclei + 2.0 *
                                    ab[IDX_HDI] / 3.0 / Hnuclei;
    IJth(A, IDX_ELEM_D, IDX_ELEM_H) = 0.0 + 1.0 * ab[IDX_HDI] / 3.0 / Hnuclei;
    IJth(A, IDX_ELEM_H, IDX_ELEM_He) = 0.0;
    IJth(A, IDX_ELEM_H, IDX_ELEM_D) = 0.0 + 2.0 * ab[IDX_HDI] / 3.0 / Hnuclei;
    IJth(A, IDX_ELEM_H, IDX_ELEM_H) = 0.0 + 4.0 * ab[IDX_H2II] / 2.0 / Hnuclei +
                                    1.0 * ab[IDX_HM] / 1.0 / Hnuclei + 1.0 *
                                    ab[IDX_HDI] / 3.0 / Hnuclei + 1.0 *
                                    ab[IDX_HI] / 1.0 / Hnuclei + 1.0 *
                                    ab[IDX_HII] / 1.0 / Hnuclei + 4.0 *
                                    ab[IDX_H2I] / 2.0 / Hnuclei;
        // clang-format on

    return NAUNET_SUCCESS;
}

// clang-format off
int RenormAbundance(realtype *rptr, realtype *ab) {
    
    ab[IDX_HeI] = ab[IDX_HeI] * (4.0 * rptr[IDX_ELEM_He] / 4.0);
    ab[IDX_HeIII] = ab[IDX_HeIII] * (4.0 * rptr[IDX_ELEM_He] / 4.0);
    ab[IDX_HeII] = ab[IDX_HeII] * (4.0 * rptr[IDX_ELEM_He] / 4.0);
    ab[IDX_H2II] = ab[IDX_H2II] * (2.0 * rptr[IDX_ELEM_H] / 2.0);
    ab[IDX_DII] = ab[IDX_DII] * (2.0 * rptr[IDX_ELEM_D] / 2.0);
    ab[IDX_DI] = ab[IDX_DI] * (2.0 * rptr[IDX_ELEM_D] / 2.0);
    ab[IDX_HM] = ab[IDX_HM] * (1.0 * rptr[IDX_ELEM_H] / 1.0);
    ab[IDX_HDI] = ab[IDX_HDI] * (2.0 * rptr[IDX_ELEM_D] / 3.0 + 1.0 * rptr[IDX_ELEM_H] / 3.0);
    ab[IDX_HI] = ab[IDX_HI] * (1.0 * rptr[IDX_ELEM_H] / 1.0);
    ab[IDX_HII] = ab[IDX_HII] * (1.0 * rptr[IDX_ELEM_H] / 1.0);
    ab[IDX_H2I] = ab[IDX_H2I] * (2.0 * rptr[IDX_ELEM_H] / 2.0);
    ab[IDX_eM] = ab[IDX_eM] * (1.0);
        // clang-format on

    return NAUNET_SUCCESS;
}