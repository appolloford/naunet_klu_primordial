// 
#include <math.h>
#include <stdio.h>

#include <algorithm>

#include "naunet_constants.h"
#include "naunet_macros.h"
#include "naunet_physics.h"
#include "naunet_utilities.h"

// clang-format off
double GetElementAbund(double *y, int elemidx) {
    if (elemidx == IDX_ELEM_He) {
        return 1.0*y[IDX_HeI] + 1.0*y[IDX_HeIII] + 1.0*y[IDX_HeII] + 0.0;
    }
    if (elemidx == IDX_ELEM_D) {
        return 1.0*y[IDX_DII] + 1.0*y[IDX_DI] + 1.0*y[IDX_HDI] + 0.0;
    }
    if (elemidx == IDX_ELEM_H) {
        return 2.0*y[IDX_H2II] + 1.0*y[IDX_HM] + 1.0*y[IDX_HDI] + 1.0*y[IDX_HI] + 
               1.0*y[IDX_HII] + 2.0*y[IDX_H2I] + 0.0;
    }
    
}

double GetMantleDens(double *y) {
    return  + 0.0;
}

double GetHNuclei(double *y) {
#ifdef IDX_ELEM_H
    return GetElementAbund(y, IDX_ELEM_H);
#else
    return 0.0;
#endif
}

double GetMu(double *y) {
    // TODO: exclude electron, grain?
    double mass = 4.0*y[IDX_HeI] + 4.0*y[IDX_HeIII] + 4.0*y[IDX_HeII] + 2.0*y[IDX_H2II] + 
                  2.0*y[IDX_DII] + 2.0*y[IDX_DI] + 1.0*y[IDX_HM] + 3.0*y[IDX_HDI] + 
                  1.0*y[IDX_HI] + 1.0*y[IDX_HII] + 2.0*y[IDX_H2I] + 0.0*y[IDX_eM] + 0.0;
    double num = y[IDX_HeI] + y[IDX_HeIII] + y[IDX_HeII] + y[IDX_H2II] +
                 y[IDX_DII] + y[IDX_DI] + y[IDX_HM] + y[IDX_HDI] + y[IDX_HI] +
                 y[IDX_HII] + y[IDX_H2I] + y[IDX_eM] + 0.0;

    return mass / num;
}

double GetGamma(double *y) {
    // TODO: different ways to get adiabatic index
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
        case 1:
            shielding = GetH2shieldingFGK(coldens);
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
        case 1:
            shielding = GetCOshieldingInt1(h2col, coldens);
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

    printf("WARNING!! Not Implemented! Return H2 shielding = -1.\n");

    /* */

    return shielding;
}

// Calculates the line self shielding function
// Ref: Federman et al. apj vol.227 p.466.
// Originally implemented in UCLCHEM
// clang-format off
double GetH2shieldingFGK(double coldens) {
    // clang-format on

    const double dopplerwidth       = 3.0e10;
    const double radiativewidth     = 8.0e7;
    const double oscillatorstrength = 1.0e-2;

    double shielding                = -1.0;

    double taud = 0.5 * coldens * 1.5e-2 * oscillatorstrength / dopplerwidth;

    // Calculate wing contribution of self shielding function sr
    if (taud < 0.0) taud = 0.0;

    double sr = 0.0;
    if (radiativewidth != 0.0) {
        double r = radiativewidth / (1.7724539 * dopplerwidth);
        double t = 3.02 * pow(1000.0 * r, -0.064);
        double u = pow(taud * r, 0.5) / t;
        sr       = pow((u * u + 0.78539816), -0.5) * r / t;
    }

    // Calculate doppler contribution of self shielding function sj
    double sj = 0.0;
    if (taud == 0.0) {
        sj = 1.0;
    } else if (taud < 2.0) {
        sj = exp(-0.6666667 * taud);
    } else if (taud < 10.0) {
        sj = 0.638 * pow(taud, -1.25);
    } else if (taud < 100.0) {
        sj = 0.505 * pow(taud, -1.15);
    } else {
        sj = 0.344 * pow(taud, -1.0667);
    }

    shielding = sj + sr;

    return shielding;
}

// Interpolate/Extropolate from table (must be rendered in naunet constants)
// clang-format off
double GetCOshieldingInt(double tgas, double h2col, double coldens) {
    // clang-format on
    double shielding = -1.0;

    /* */

    printf("WARNING!! Not Implemented! Return CO shielding = -1.\n");

    /* */

    return shielding;
}

// clang-format off
double GetCOshieldingInt1(double h2col, double coldens) {
    // clang-format on
    double shielding = -1.0;

    /* */

    printf("WARNING!! Not Implemented! Return CO shielding = -1.\n");

    /* */

    return shielding;
}

// Interpolate/Extropolate from table (must be rendered in naunet constants)
// clang-format off
double GetN2shieldingInt(double tgas, double h2col, double coldens) {
    // clang-format on

    double shielding = -1.0;

    /* */

    printf("WARNING!! Not Implemented! Return N2 shielding = -1.\n");

    /* */

    return shielding;
}

// Calculate xlamda := tau(lambda) / tau(visual)
// tau(lambda) is the opt. depth for dust extinction at
// wavelength x (cf. b.d.savage and j.s.mathis, annual
// review of astronomy and astrophysics vol.17(1979),p.84)
// clang-format off
double xlamda(double wavelength) {
    // clang-format on
    double x[29] = {910.0,  950.0,  1000.0,  1050.0,  1110.0, 1180.0,
                    1250.0, 1390.0, 1490.0,  1600.0,  1700.0, 1800.0,
                    1900.0, 2000.0, 2100.0,  2190.0,  2300.0, 2400.0,
                    2500.0, 2740.0, 3440.0,  4000.0,  4400.0, 5500.0,
                    7000.0, 9000.0, 12500.0, 22000.0, 34000.0};

    double y[29] = {5.76, 5.18, 4.65, 4.16, 3.73, 3.4,  3.11, 2.74, 2.63, 2.62,
                    2.54, 2.5,  2.58, 2.78, 3.01, 3.12, 2.86, 2.58, 2.35, 2.0,
                    1.58, 1.42, 1.32, 1.0,  0.75, 0.48, 0.28, 0.12, 0.05};

    if (wavelength < x[0]) {
        return 5.76;
    }

    else if (wavelength >= x[28]) {
        return 0.05 - 5.16e-11 * (wavelength - x[28]);
    }

    for (int i = 0; i < 28; i++) {
        if (wavelength >= x[i] && wavelength < x[i + 1]) {
            return y[i] +
                   (y[i + 1] - y[i]) * (wavelength - x[i]) / (x[i + 1] - x[i]);
        }
    }

    return 0.0;
}

// Calculate the influence of dust extinction (g=0.8, omega=0.3)
// Ref: Wagenblast & Hartquist, mnras237, 1019 (1989)
// Adapted from UCLCHEM
// clang-format off
double GetGrainScattering(double av, double wavelength) {
    // clang-format on
    double c[6] = {1.0e0, 2.006e0, -1.438e0, 7.364e-1, -5.076e-1, -5.920e-2};
    double k[6] = {7.514e-1, 8.490e-1, 1.013e0, 1.282e0, 2.005e0, 5.832e0};

    double tv   = av / 1.086;
    double tl   = tv * xlamda(wavelength);

    double scat = 0.0;
    double expo;
    if (tl < 1.0) {
        expo = k[0] * tl;
        if (expo < 35.0) {
            scat = c[0] * exp(-expo);
        }
    } else {
        for (int i = 1; i < 6; i++) {
            expo = k[i] * tl;
            if (expo < 35.0) {
                scat = scat + c[i] * exp(-expo);
            }
        }
    }

    return scat;
}

// Calculate lambda bar (in a) according to equ. 4 of van dishoeck
// and black, apj 334, p771 (1988)
// Adapted from UCLCHEM
// clang-format off
double GetCharactWavelength(double h2col, double cocol) {
    // clang-format on
    double logco = log10(abs(cocol) + 1.0);
    double logh2 = log10(abs(h2col) + 1.0);

    double lbar  = (5675.0 - 200.6 * logh2) - (571.6 - 24.09 * logh2) * logco +
                  (18.22 - 0.7664 * logh2) * pow(logco, 2.0);

    // lbar represents the mean of the wavelengths of the 33
    // dissociating bands weighted by their fractional contribution
    // to the total rate of each depth. lbar cannot be larger than
    // the wavelength of band 33 (1076.1a) and not be smaller than
    // the wavelength of band 1 (913.6a).

    /* */
    lbar = std::min(1076.0, std::max(913.0, lbar));
    /* */
    return lbar;
}