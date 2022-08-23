// 
#include <stdio.h>

#include <stdexcept>
#include <vector>

#include "naunet.h"
#include "naunet_data.h"
#include "naunet_macros.h"
#include "naunet_ode.h"
#include "naunet_timer.h"

int main() {
    double spy = 86400.0 * 365.0;

    NaunetData data;
    //
    double nH       = 1e-2;
    double Tgas     = 100000.0;

    data.nH         = nH;
    data.Tgas       = Tgas;


    Naunet naunet;
    if (naunet.Init() == NAUNET_FAIL) {
        printf("Initialize Fail\n");
        return 1;
    }

#ifdef USE_CUDA
    if (naunet.Reset(1) == NAUNET_FAIL) {
        throw std::runtime_error("Fail to reset the number of systems");
    }
#endif

    //
    double y[NEQUATIONS];
    for (int i = 0; i < NEQUATIONS; i++) {
        y[i] = 1.e-40;
    }
    y[IDX_HI]   = nH;
    y[IDX_HII]  = 1e-4 * nH;
    y[IDX_HeI]  = 1e-1 * nH;
    y[IDX_HDI]  = 1.5e-5 * nH;
    y[IDX_H2I]  = 1.5e-5 * nH;
    y[IDX_eM]   = 1e-4 * nH;
    y[IDX_TGAS] = Tgas;


    FILE *fbin = fopen("evolution_singlegrid.bin", "w");
    FILE *ftxt = fopen("evolution_singlegrid.txt", "w");
    FILE *ttxt = fopen("time_singlegrid.txt", "w");
#ifdef NAUNET_DEBUG
    printf("Initialization is done. Start to evolve.\n");
    FILE *rtxt               = fopen("reactionrates.txt", "w");
    double rates[NREACTIONS] = {0.0};
#endif

    //
    std::vector<double> timesteps;
    double logtstart = 3.0, logtend = 7.0, logtstep = 0.1;
    double time = 0.0;
    for (double logtime = logtstart; logtime < logtend + 0.1 * logtstep;
         logtime += logtstep) {
        double dtyr = pow(10.0, logtime) - time;
        timesteps.push_back(dtyr);
        time += dtyr;
    }
    //

    double dtyr = 0.0, curtime = 0.0;

    // write the initial abundances
    fwrite(&curtime, sizeof(double), 1, fbin);
    fwrite(y, sizeof(double), NEQUATIONS, fbin);

    fprintf(ftxt, "%13.7e ", curtime);
    for (int j = 0; j < NEQUATIONS; j++) {
        fprintf(ftxt, "%13.7e ", y[j]);
    }
    fprintf(ftxt, "\n");

    for (auto step = timesteps.begin(); step != timesteps.end(); step++) {
#ifdef NAUNET_DEBUG
        EvalRates(rates, y, &data);
        for (int j = 0; j < NREACTIONS; j++) {
            fprintf(rtxt, "%13.7e ", rates[j]);
        }
        fprintf(rtxt, "\n");
#endif

        //
        // synchronize the temperature if you want to make chemistry and
        // heating/cooling consistent
        data.Tgas = y[IDX_TGAS];


        dtyr = *step;

        Timer timer;
        timer.start();
        naunet.Solve(y, dtyr * spy, &data);
        timer.stop();

        curtime += dtyr;

        // write the abundances after each step
        fwrite(&curtime, sizeof(double), 1, fbin);
        fwrite(y, sizeof(double), NEQUATIONS, fbin);

        fprintf(ftxt, "%13.7e ", curtime);
        for (int j = 0; j < NEQUATIONS; j++) {
            fprintf(ftxt, "%13.7e ", y[j]);
        }
        fprintf(ftxt, "\n");

        // float duration = (float)timer.elapsed() / 1e6;
        double duration = timer.elapsed();
        fprintf(ttxt, "%8.5e \n", duration);
        printf("Time = %13.7e yr, elapsed: %8.5e sec\n", curtime, duration);
    }

    fclose(fbin);
    fclose(ftxt);
    fclose(ttxt);
#ifdef NAUNET_DEBUG
    fclose(rtxt);
#endif

    if (naunet.Finalize() == NAUNET_FAIL) {
        printf("Finalize Fail\n");
        return 1;
    }

    return 0;
}