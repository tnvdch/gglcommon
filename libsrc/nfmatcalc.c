/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/nfmatcalc.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations on 3x3 matrix and 3x1 vectors
 *            (regroups routines taken from numerical recipies
 *            and adopted to struct matrix as defined in
 *            matcalc.h)
 *
 * $Id: nfmatcalc.c,v 1.6 1999/09/08 09:12:18 favey Exp $
 *
 *    NOTICE: This is commercial software and may NOT be distributed
 *            nor be used outside of the Geodesy and Geodynamics Lab, ETHZ
 *
 *   AUTHORS: Cocard Marc     (adopted from numerical recipies)
 *
 *   HISTORY: 
 * $Log: nfmatcalc.c,v $
 * Revision 1.6  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.5  1998/11/03 18:00:51  favey
 * hardcoded double instead of MATTYPE, as changes go now into sed_matcalc
 *
 * Revision 1.4  1998/06/24 09:41:50  cocard
 * *** empty log message ***
 *
 * Revision 1.3  1998/06/23 17:12:51  cocard
 * *** empty log message ***
 *
 * Revision 1.2  1998/06/23 17:06:40  cocard
 * *** empty log message ***
 *
 *
 */

/*
 *  Standard Libraries
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 *  our special headers
 */
#include "matcalc.h"
#include "nfmatcalc.h"

/* routine m_jacobi taken from
 *
 * (C) Copr. 1986-92 Numerical Recipes Software '',.
 *
 * adopted to struct matrix
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

int m_jacobi1(struct matrix *a, struct matrix *eval, struct matrix *evec, int *nrot)
{

#define MAXITER 500
#define ROTATE(a,i,j,k,l) g=a->mat[i][j];h=a->mat[k][l];\
                          a->mat[i][j]=g-s*(h+g*tau);\
                          a->mat[k][l]=h+s*(g-h*tau);

  int i,j,iq,ip;
  int iter, n;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
  
  n=a->row;
  b = (double *) malloc(n*sizeof(double));
  z = (double *) malloc(n*sizeof(double));

#ifdef WITH_DEBUG
  if (b==NULL || z==NULL)
    {
      m_error("m_jacobi",NULL, NULL, NULL, NULL, NULL);
      fprintf(stderr,"Numerical Recipes run-time error...\n");
      fprintf(stderr,"malloc error\n");
      return -1; 
    }
#endif

  m_unit(evec);

  for (i=0;i<n;i++) { 
   b[i]=eval->mat[0][i]=a->mat[i][i];
    z[i]=0.0;
  }

  *nrot=0;

  for (iter=1;iter<=MAXITER; iter++) {

    sm=0.0;

    for (ip=0; ip<n-1; ip++) {  
      for (iq=ip+1;iq<n; iq++)
	sm += fabs(a->mat[ip][iq]);
    }

    if (sm == 0.0)  /* successful computation */
      { 
	/* refill the upper part of a matrix symmetrically */
	for (i=0; i<n-1; i++)
	  for (j=i; j<n; j++)
	    a->mat[i][j]=a->mat[j][i];
	free(z);
	free(b);
	return 0;
      }

    if (iter < 4) 
      tresh=0.2*sm/(n*n);
    else 
      tresh=0.0;

    for (ip=0;ip<n-1;ip++)      
      for (iq=ip+1;iq<n;iq++)   
	{
	  g=100.0*fabs(a->mat[ip][iq]);

	  if (iter > 4 && 
	      (double)(fabs(eval->mat[0][ip])+g) == (double)fabs(eval->mat[0][ip]) && 
	      (double)(fabs(eval->mat[0][iq])+g) == (double)fabs(eval->mat[0][iq]))
	    a->mat[ip][iq]=0.0;
	  else if (fabs(a->mat[ip][iq]) > tresh) {
	    h=eval->mat[0][iq] - eval->mat[0][ip];

	    if ((double)(fabs(h)+g) == (double)fabs(h))
	      t=a->mat[ip][iq]/h;
	    else {
	      theta=0.5*h/(a->mat[ip][iq]);
	      t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	      if (theta < 0.0) t = -t;
	    }
	    
	    c=1.0/sqrt(1+t*t);
	    s=t*c;
	    tau=s/(1.0+c);
	    h=t*a->mat[ip][iq];
	    z[ip] -= h;
	    z[iq] += h;
	    eval->mat[0][ip] -= h;
	    eval->mat[0][iq] += h;
	    a->mat[ip][iq]=0.0;

	    for (j=0;    j<ip; j++) { ROTATE(a,j,ip,j,iq) };
	    for (j=ip+1; j<iq; j++) { ROTATE(a,ip,j,j,iq) };	   
	    for (j=iq+1; j<n;  j++) { ROTATE(a,ip,j,iq,j) };
	    for (j=0; j<n; j++) { ROTATE(evec,j,ip,j,iq) };

	    ++(*nrot);
	  }	  
	}

    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      eval->mat[0][ip]=b[ip];
      z[ip]=0.0;
    }
  }
  free(z);
  free(b);
    
#ifdef WITH_DEBUG
  m_error("m_jacobi",NULL, NULL, NULL, NULL, NULL);
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n","Too many iterations in routine jacobi");
#endif
 
  return -1;
}


/*
 * difference to m_jacobi1 : matrix a is replaced by its eigenvalue matrix.
 */
int m_jacobi2(struct matrix *a, struct matrix *evec, int *nrot)
{
  struct matrix *eval;
  int err, i;

  eval=m_init(1, a->row);
  if (eval==NULL) return -1;
  err=m_jacobi1(a, eval, evec, nrot);
  m_zero(a);
  for (i=0; i<a->row; i++) a->mat[i][i]=eval->mat[0][i];
  m_free(eval);
  return err;
}

