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
    int nsystem      = 2048;
    double spy       = 86400.0 * 365.0;

    NaunetData *data = new NaunetData[nsystem];
    //
    double nH       = 1e-2;
    double Tgas     = 100000.0;

    for (int isys = 0; isys < nsystem; isys++) {
        data[isys].nH       = nH;
        data[isys].Tgas     = Tgas;
    }


    Naunet naunet;
    if (naunet.Init() == NAUNET_FAIL) {
        printf("Initialize Fail\n");
        return 1;
    }

    if (naunet.Reset(nsystem) == NAUNET_FAIL) {
        throw std::runtime_error("Fail to reset the number of systems");
    }

    //
    double *y = new double[nsystem * NEQUATIONS];
    for (int isys = 0; isys < nsystem; isys++) {
        for (int i = 0; i < NEQUATIONS; i++) {
            y[isys * NEQUATIONS + i] = 1.e-40;
        }
        y[isys * NEQUATIONS + IDX_HI]   = nH;
        y[isys * NEQUATIONS + IDX_HII]  = 1e-4 * nH;
        y[isys * NEQUATIONS + IDX_HeI]  = 1e-1 * nH;
        y[isys * NEQUATIONS + IDX_HDI]  = 1.5e-5 * nH;
        y[isys * NEQUATIONS + IDX_H2I]  = 1.5e-5 * nH;
        y[isys * NEQUATIONS + IDX_eM]   = 1e-4 * nH;
        y[isys * NEQUATIONS + IDX_TGAS] = Tgas;
    }


    FILE *fbin = fopen("evolution_multiplegrid.bin", "w");
    FILE *ftxt = fopen("evolution_multiplegrid.txt", "w");
    FILE *ttxt = fopen("time_parallel.txt", "w");

#ifdef NAUNET_DEBUG
    printf("Initialization is done. Start to evolve.\n");
    // FILE *rtxt = fopen("reactionrates.txt", "w");
    // double rates[NREACTIONS];
#endif

    //
    std::vector<double> timesteps;
    double logtstart = 2.0, logtend = 4.0, logtstep = 0.1;
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
    for (int isys = 0; isys < nsystem; isys++) {
        fwrite((double *)&isys, sizeof(double), 1, fbin);
        fwrite(&curtime, sizeof(double), 1, fbin);
        fwrite(&y[isys * NEQUATIONS], sizeof(double), NEQUATIONS, fbin);

        fprintf(ftxt, "%13.7e ", (double)isys);
        fprintf(ftxt, "%13.7e ", curtime);
        for (int j = 0; j < NEQUATIONS; j++) {
            fprintf(ftxt, "%13.7e ", y[isys * NEQUATIONS + j]);
        }
        fprintf(ftxt, "\n");
    }
    for (auto step = timesteps.begin(); step != timesteps.end(); step++) {
#ifdef NAUNET_DEBUG
        // EvalRates only receive one system as input, disabled in parallel test
        // EvalRates(rates, y, data);
        // for (int j = 0; j < NREACTIONS; j++) {
        //     fprintf(rtxt, "%13.7e ", rates[j]);
        // }
        // fprintf(rtxt, "\n");
#endif

        //
        // synchronize the temperature if you want to make chemistry and
        // heating/cooling consistent
        for (int isys = 0; isys < nsystem; isys++) {
            data[isys].Tgas = y[isys * NEQUATIONS + IDX_TGAS];
        }


        dtyr = *step;

        Timer timer;
        timer.start();
        naunet.Solve(y, dtyr * spy, data);
        timer.stop();

        curtime += dtyr;

        // write the abundances after each step
        for (int isys = 0; isys < nsystem; isys++) {
            fwrite((double *)&isys, sizeof(double), 1, fbin);
            fwrite(&curtime, sizeof(double), 1, fbin);
            fwrite(&y[isys * NEQUATIONS], sizeof(double), NEQUATIONS, fbin);

            fprintf(ftxt, "%13.7e ", (double)isys);
            fprintf(ftxt, "%13.7e ", curtime);
            for (int j = 0; j < NEQUATIONS; j++) {
                fprintf(ftxt, "%13.7e ", y[isys * NEQUATIONS + j]);
            }
            fprintf(ftxt, "\n");
        }

        // float duration = (float)timer.elapsed() / 1e6;
        double duration = timer.elapsed();
        fprintf(ttxt, "%8.5e \n", duration);
        printf("Time = %13.7e yr, elapsed: %8.5e sec\n", curtime, duration);
    }

    fclose(fbin);
    fclose(ftxt);
    fclose(ttxt);

#ifdef NAUNET_DEBUG
    // fclose(rtxt);
#endif

    if (naunet.Finalize() == NAUNET_FAIL) {
        printf("Finalize Fail\n");
        return 1;
    }

    delete[] data;
    delete[] y;

    return 0;
}