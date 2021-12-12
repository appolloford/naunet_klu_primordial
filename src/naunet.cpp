#include <cvode/cvode.h>  // prototypes for CVODE fcts., consts.
/* */
#include <nvector/nvector_serial.h>      // access to serial N_Vector
#include <sunlinsol/sunlinsol_klu.h>     // access to KLU sparse direct solver
#include <sunmatrix/sunmatrix_sparse.h>  // access to sparse SUNMatrix
/* */
/*  */
#include "naunet.h"
/*  */
#include "naunet_ode.h"

// check_flag function is from the cvDiurnals_ky.c example from the CVODE
// package. Check function return value...
//   opt == 0 means SUNDIALS function allocates memory so check if
//            returned NULL pointer
//   opt == 1 means SUNDIALS function returns a flag so check if
//            flag >= 0
//   opt == 2 means function allocates memory so check if returned
//            NULL pointer
static int check_flag(void *flagvalue, const char *funcname, int opt) {
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
        fprintf(stderr,
                "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return 1;
    }

    /* Check if flag < 0 */
    else if (opt == 1) {
        errflag = (int *)flagvalue;
        if (*errflag < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                    funcname, *errflag);
            return 1;
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL) {
        fprintf(stderr,
                "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return 1;
    }

    return 0;
}

Naunet::Naunet(){};

Naunet::~Naunet(){};

int Naunet::Init(int nsystem, double atol, double rtol) {
    n_system_ = nsystem;
    atol_     = atol;
    rtol_     = rtol;

    /* */
    if (nsystem != 1) {
        printf("This solver doesn't support nsystem > 1!");
        return NAUNET_FAIL;
    }

    cv_y_  = N_VNew_Serial((sunindextype)NEQUATIONS);
    cv_a_  = SUNSparseMatrix(NEQUATIONS, NEQUATIONS, NNZ, CSR_MAT);
    cv_ls_ = SUNLinSol_KLU(cv_y_, cv_a_);

    /*  */

    cv_mem_ = CVodeCreate(CV_BDF);

    int flag;
    flag = CVodeInit(cv_mem_, Fex, 0.0, cv_y_);
    if (check_flag(&flag, "CVodeInit", 1)) return 1;
    flag = CVodeSStolerances(cv_mem_, rtol_, atol_);
    if (check_flag(&flag, "CVodeSStolerances", 1)) return 1;
    flag = CVodeSetLinearSolver(cv_mem_, cv_ls_, cv_a_);
    if (check_flag(&flag, "CVodeSetLinearSolver", 1)) return 1;
    flag = CVodeSetJacFn(cv_mem_, Jac);
    if (check_flag(&flag, "CVodeSetJacFn", 1)) return 1;

    // reset the n_vector to empty, maybe not necessary
    /* */

    // N_VDestroy(cv_y_);
    // cv_y_ = N_VNewEmpty_Serial((sunindextype)NEQUATIONS);

    /* */

    return NAUNET_SUCCESS;
};

int Naunet::Finalize() {
    // N_VDestroy(cv_y_);
    N_VFreeEmpty(cv_y_);
    SUNMatDestroy(cv_a_);
    CVodeFree(&cv_mem_);
    SUNLinSolFree(cv_ls_);
    // delete m_data;

    /*  */

    return NAUNET_SUCCESS;
};

/*  */

int Naunet::Solve(realtype *ab, realtype dt, NaunetData *data) {
    realtype t0 = 0.0;
    int flag;

    /* */

    // realtype *ydata = N_VGetArrayPointer(cv_y_);
    // for (int i=0; i<NEQUATIONS; i++)
    // {
    //     ydata[i] = ab[i];
    // }
    N_VSetArrayPointer(ab, cv_y_);

    /*  */

#ifdef NAUNET_DEBUG
    /*  */
#endif

    flag = CVodeReInit(cv_mem_, 0.0, cv_y_);
    if (check_flag(&flag, "CVodeReInit", 1)) return 1;
    flag = CVodeSetUserData(cv_mem_, data);
    if (check_flag(&flag, "CVodeSetUserData", 1)) return 1;

    flag = CVode(cv_mem_, dt, cv_y_, &t0, CV_NORMAL);

    /* */

    ab   = N_VGetArrayPointer(cv_y_);

    /* */

    return NAUNET_SUCCESS;
};

#ifdef PYMODULE
py::array_t<realtype> Naunet::PyWrapSolve(py::array_t<realtype> arr,
                                          realtype dt, NaunetData *data) {
    py::buffer_info info = arr.request();
    realtype *ab         = static_cast<realtype *>(info.ptr);

    Solve(ab, dt, data);

    return py::array_t<realtype>(info.shape, ab);
}
#endif