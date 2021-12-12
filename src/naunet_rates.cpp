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

// clang-format off
int EvalRates(realtype *k, realtype *y, NaunetData *u_data) {

    realtype nH = u_data->nH;
    realtype Tgas = u_data->Tgas;
    
    
    // clang-format on

    // Some variable definitions from krome
    realtype Te      = Tgas * 8.617343e-5;            // Tgas in eV (eV)
    realtype lnTe    = log(Te);                       // ln of Te (#)
    realtype T32     = Tgas * 0.0033333333333333335;  // Tgas/(300 K) (#)
    realtype invT    = 1.0 / Tgas;                    // inverse of T (1/K)
    realtype invTe   = 1.0 / Te;                      // inverse of T (1/eV)
    realtype sqrTgas = sqrt(Tgas);  // Tgas rootsquare (K**0.5)

    // reaaction rate (k) of each reaction
    // clang-format off
    k[0] = exp(-32.71396786e0 + 13.5365560e0*lnTe-5.73932875e0*(pow(lnTe,
        2)) + 1.56315498e0*(pow(lnTe, 3)) - 0.28770560e0*(pow(lnTe, 4)) +
        3.48255977e-2*(pow(lnTe, 5)) - 2.63197617e-3*(pow(lnTe, 6)) +
        1.11954395e-4*(pow(lnTe, 7)) - 2.03914985e-6*(pow(lnTe, 8)));
        
    if (Tgas>1.0 && Tgas<5500.0) { k[1] = 3.92e-13*pow(invTe, 0.6353e0);  }
        
    k[2] = exp(-28.61303380689232e0 -
        0.7241125657826851e0*lnTe-0.02026044731984691e0*pow(lnTe, 2) -
        0.002380861877349834e0*pow(lnTe, 3) - 0.0003212605213188796e0*pow(lnTe,
        4) - 0.00001421502914054107e0*pow(lnTe, 5) +
        4.989108920299513e-6*pow(lnTe, 6) + 5.755614137575758e-7*pow(lnTe, 7) -
        1.856767039775261e-8*pow(lnTe, 8) - 3.071135243196595e-9*pow(lnTe, 9));
        
    k[3] = exp(-44.09864886e0 + 23.91596563e0*lnTe-10.7532302e0*(pow(lnTe,
        2)) + 3.05803875e0*(pow(lnTe, 3)) - 0.56851189e0*(pow(lnTe, 4)) +
        6.79539123e-2*(pow(lnTe, 5)) - 5.00905610e-3*(pow(lnTe, 6)) +
        2.06723616e-4*(pow(lnTe, 7)) - 3.64916141e-6*(pow(lnTe, 8)));
        
    if (Tgas>1.0 && Tgas<9280.0) { k[4] = 3.92e-13*pow(invTe, 0.6353e0);  }
        
    k[5] = 1.54e-9*(1.e0 +
        0.3e0/exp(8.099328789667e0*invTe))/(exp(40.49664394833662e0*invTe)*pow(Te,
        1.5e0)) + 3.92e-13/pow(Te, 0.6353e0);
        
    k[6] = exp(-68.71040990212001e0 +
        43.93347632635e0*lnTe-18.48066993568e0*pow(lnTe, 2) +
        4.701626486759002e0*pow(lnTe, 3) - 0.7692466334492e0*pow(lnTe, 4) +
        0.08113042097303e0*pow(lnTe, 5) - 0.005324020628287001e0*pow(lnTe, 6) +
        0.0001975705312221e0*pow(lnTe, 7) - 3.165581065665e-6*pow(lnTe, 8));
        
    k[7] = 3.36e-10/sqrTgas/pow((Tgas/1.e3), 0.2e0)/(1 + pow((Tgas/1.e6),
        0.7e0));
        
    k[8] = 6.77e-15*pow(Te, 0.8779e0);
        
    if (Tgas>1.0 && Tgas<1160.0) { k[9] = 1.43e-9;  }
        
    k[10] = exp(-20.06913897587003e0 +
        0.2289800603272916e0*lnTe+0.03599837721023835e0*pow(lnTe, 2) -
        0.004555120027032095e0*pow(lnTe, 3) - 0.0003105115447124016e0*pow(lnTe,
        4) + 0.0001073294010367247e0*pow(lnTe, 5) -
        8.36671960467864e-6*pow(lnTe, 6) + 2.238306228891639e-7*pow(lnTe, 7));
        
    if (Tgas>1.0 && Tgas<6700.0) { k[11] = 1.85e-23*pow(Tgas, 1.8e0);  }
        
    k[12] = 5.81e-16*pow((Tgas/5.62e4), (-0.6657e0*log10(Tgas/5.62e4)));
        
    k[13] = 6.0e-10;
        
    k[14] = exp(-24.24914687731536e0 +
        3.400824447095291e0*lnTe-3.898003964650152e0*pow(lnTe, 2) +
        2.045587822403071e0*pow(lnTe, 3) - 0.5416182856220388e0*pow(lnTe, 4) +
        0.0841077503763412e0*pow(lnTe, 5) - 0.007879026154483455e0*pow(lnTe, 6)
        + 0.0004138398421504563e0*pow(lnTe, 7) - 9.36345888928611e-6*pow(lnTe,
        8));
        
    k[15] = 5.6e-11*exp(-102124.e0*invT)*pow(Tgas, 0.5e0);
        
    k[16] = 1.0670825e-10*pow(Te, 2.012e0)*exp(-4.463e0*invTe)/pow((1.e0 +
        0.2472e0*Te), 3.512e0);
        
    k[17] = exp(-18.01849334273e0 +
        2.360852208681e0*lnTe-0.2827443061704e0*pow(lnTe, 2) +
        0.01623316639567e0*pow(lnTe, 3) - 0.03365012031362999e0*pow(lnTe, 4) +
        0.01178329782711e0*pow(lnTe, 5) - 0.001656194699504e0*pow(lnTe, 6) +
        0.0001068275202678e0*pow(lnTe, 7) - 2.631285809207e-6*pow(lnTe, 8));
        
    if (Tgas>1.0 && Tgas<1160.0) { k[18] = 2.56e-9*pow(Te, 1.78186e0);  }
        
    k[19] = exp(-20.37260896533324e0 +
        1.139449335841631e0*lnTe-0.1421013521554148e0*pow(lnTe, 2) +
        0.00846445538663e0*pow(lnTe, 3) - 0.0014327641212992e0*pow(lnTe, 4) +
        0.0002012250284791e0*pow(lnTe, 5) + 0.0000866396324309e0*pow(lnTe, 6) -
        0.00002585009680264e0*pow(lnTe, 7) + 2.4555011970392e-6*pow(lnTe, 8) -
        8.06838246118e-8*pow(lnTe, 9));
        
    k[20] = 6.5e-9/sqrt(Te);
        
    k[21] = 1.e-8*pow(Tgas, (-0.4e0));
        
    if (Tgas>1.0 && Tgas<617.0) { k[22] = 1.e-8;  }
        
    k[23] = 1.32e-6*pow(Tgas, (-0.76e0));
        
    k[24] = 5.e-7*sqrt(1.e2*invT);
        
    if (Tgas>1.0 && Tgas<300.0) { k[25] = 1.3e-32*pow((T32), (-0.38e0));  }
        
    k[26] = 1.3e-32*pow((T32), (-1.00e0));
        
    if (Tgas>1.0 && Tgas<300.0) { k[27] = 1.3e-32*pow((T32),
        (-0.38e0))/8.e0;  }
        
    k[28] = 1.3e-32*pow((T32), (-1.00e0))/8.e0;
        
    k[29] = 2.00e-10*pow(Tgas, (0.402e0))*exp(-37.1e0*invT) -
        3.31e-17*pow(Tgas, (1.48e0));
        
    k[30] = 2.06e-10*pow(Tgas, (0.396))*exp(-33.e0*invT) + 2.03e-9*pow(Tgas,
        (-0.332));
        
    k[31] = 1.e-9*(0.417 + 0.846*log10(Tgas) - 0.137*pow((log10(Tgas)), 2));
        
    k[32] = 1.0e-9*exp(-4.57e2*invT);
        
    if (Tgas>1.0 && Tgas<2000.0) { k[33] = pow(10, (-56.4737 +
        5.88886*log10(Tgas) + 7.19692*pow((log10(Tgas)), 2) +
        2.25069*pow((log10(Tgas)), 3) - 2.16903*pow((log10(Tgas)), 4) +
        0.317887*pow((log10(Tgas)), 5)));  }
        
    k[34] = 3.17e-10*exp(-5207.*invT);
        
    k[35] = 5.25e-11*exp(-4430.*invT + 1.739e5*pow((invT), 2));
        
    k[36] = 1.5e-9*pow((T32), (-0.1e0));
        
    k[37] = 3.6e-12*pow((Tgas/300), (-0.75e0));
        
    
        // clang-format on

    return NAUNET_SUCCESS;
}