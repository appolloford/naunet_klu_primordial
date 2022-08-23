// 
#ifndef __NAUNET_DATA_H__
#define __NAUNET_DATA_H__

// Struct for holding the nessesary additional variables for the problem.
struct NaunetData {
    // clang-format off
    double nH;
    double Tgas;
    double mu = -1.000e+00;
    double gamma = -1.000e+00;
    
    // clang-format on
};

#endif