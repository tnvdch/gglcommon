/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/nfmatcalc.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations (non free stuff taken from numerical recipies)
 *
 * $Id: nfmatcalc.h,v 1.3 1998/11/02 10:42:44 favey Exp $
 *
 *   AUTHORS: Cocard Marc 
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            23-06-1998  --  created
 * $Log: nfmatcalc.h,v $
 * Revision 1.3  1998/11/02 10:42:44  favey
 * added RCS verbosity
 *
 *
 *   CAVEATS: 
 *
 */

#ifndef NF_MATCALC
#define NF_MATCALC
#ifndef GGL_MATCALC
*** ERROR: you must include matcalc.h first, nfmatcalc relies on it!
#endif
/*
 *  JACOBI : calculates Eigenvalues and Eigenvectors from a symmetric matrix a
 *
 * m_jacobi1 : input  : matrix a is conserved
 *             output : eigenvalues (eval) stored in a vector (matrix 1 x DIM)
 *             output : evectors (evec) in form of a rotation matrix
 *             output : total number of individual rotations during  iteration
 *             return value = error
 *
 * m_jacobi2 : input/output : matrix a is overwritten. New matrix a is diagonal
 *                            and contains the eigenvalues (A replaced by D so 
 *                            that a = ev * d * evT)
 *             output : evectors (evec) in form of a rotation matrix
 *             output : total number of individual rotations during  iteration
 *             return value = error
 */

int m_jacobi1(struct matrix *a, struct matrix *eval, struct matrix *evec, int *nrot);

int m_jacobi2(struct matrix *a, struct matrix *evec, int *nrot);

#endif  /* NF_MATCALC */   
