#include <math.h>

#include "naunet_constants.h"
#include "naunet_macros.h"
#include "naunet_physics.h"

// clang-format off
double GetMantleDens(double *y) {
    return 0.0;
}

double GetMu(double *y) {
    return (y[IDX_DI]*2.0 + y[IDX_DII]*2.0 + y[IDX_GRAINI]*0.0 + y[IDX_HI]*1.0 +
        y[IDX_HII]*1.0 + y[IDX_HM]*1.0 + y[IDX_H2I]*2.0 + y[IDX_H2II]*2.0 +
        y[IDX_HDI]*3.0 + y[IDX_HeI]*4.0 + y[IDX_HeII]*4.0 + y[IDX_HeIII]*4.0 +
        y[IDX_eM]*0.0) / (y[IDX_DI] + y[IDX_DII] + y[IDX_GRAINI] + y[IDX_HI] +
        y[IDX_HII] + y[IDX_HM] + y[IDX_H2I] + y[IDX_H2II] + y[IDX_HDI] +
        y[IDX_HeI] + y[IDX_HeII] + y[IDX_HeIII] + y[IDX_eM]);
}

double GetGamma(double *y) {
    return 5.0 / 3.0;
}

double GetNumDens(double *y) {
    double numdens = 0.0;

    for (int i = 0; i < NSPECIES; i++) numdens += y[i];
    return numdens;
}
// clang-format on

// clang-format off
double GetShieldingFactor(int specidx, double h2coldens, double spcoldens,
                          double tgas, int method) {
    // clang-format on
    double factor;
#ifdef IDX_H2I
    if (specidx == IDX_H2I) {
        factor = GetH2shielding(h2coldens, method);
    }
#endif
#ifdef IDX_COI
    if (specidx == IDX_COI) {
        factor = GetCOshielding(tgas, h2coldens, spcoldens, method);
    }
#endif
#ifdef IDX_N2I
    if (specidx == IDX_N2I) {
        factor = GetN2shielding(tgas, h2coldens, spcoldens, method);
    }
#endif

    return factor;
}

// clang-format off
double GetH2shielding(double coldens, int method) {
    // clang-format on
    double shielding = -1.0;
    switch (method) {
        case 0:
            shielding = GetH2shieldingInt(coldens);
            break;
        default:
            break;
    }
    return shielding;
}

// clang-format off
double GetCOshielding(double tgas, double h2col, double coldens, int method) {
    // clang-format on
    double shielding = -1.0;
    switch (method) {
        case 0:
            shielding = GetCOshieldingInt(tgas, h2col, coldens);
            break;
        default:
            break;
    }
    return shielding;
}

// clang-format off
double GetN2shielding(double tgas, double h2col, double coldens, int method) {
    // clang-format on
    double shielding = -1.0;
    switch (method) {
        case 0:
            shielding = GetN2shieldingInt(tgas, h2col, coldens);
            break;
        default:
            break;
    }
    return shielding;
}

// Interpolate/Extropolate from table (must be rendered in naunet constants)
// clang-format off
double GetH2shieldingInt(double coldens) {
    // clang-format on

    double shielding = -1.0;

    /* */

    return shielding;
}

// Interpolate/Extropolate from table (must be rendered in naunet constants)
// clang-format off
double GetCOshieldingInt(double tgas, double h2col, double coldens) {
    // clang-format on
    double shielding = -1.0;

    /* */

    return shielding;
}

// Interpolate/Extropolate from table (must be rendered in naunet constants)
// clang-format off
double GetN2shieldingInt(double tgas, double h2col, double coldens) {
    // clang-format on

    double shielding = -1.0;

    /* */

    return shielding;
}