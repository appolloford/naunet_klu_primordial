// 
#ifndef __NAUNET_MACROS_H__
#define __NAUNET_MACROS_H__

// clang-format off
#define NAUNET_SUCCESS 0
#define NAUNET_FAIL 1

#define MAX_NSYSTEMS 1

#define NELEMENTS 3
#define NSPECIES 12
#define NHEATPROCS 0
#define NCOOLPROCS 11
#define THERMAL (NHEATPROCS || NCOOLPROCS)
#if (NSPECIES + THERMAL)
#define NEQUATIONS (NSPECIES + THERMAL)
#else
#define NEQUATIONS 1
#endif
#define NREACTIONS 38
// non-zero terms in jacobian matrix, used in sparse matrix
#define NNZ 

#define IDX_ELEM_D 0
#define IDX_ELEM_H 1
#define IDX_ELEM_He 2

#define IDX_DI 0
#define IDX_DII 1
#define IDX_HI 2
#define IDX_HII 3
#define IDX_HM 4
#define IDX_H2I 5
#define IDX_H2II 6
#define IDX_HDI 7
#define IDX_HeI 8
#define IDX_HeII 9
#define IDX_HeIII 10
#define IDX_eM 11

#if THERMAL
#define IDX_TGAS NSPECIES
#endif

#endif