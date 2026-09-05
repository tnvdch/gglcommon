/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/tools.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: 
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: tools.c,v 1.7 2001/08/31 08:28:41 cocard Exp $
 *
 *  DISCLAIMER:   This program is free software; you can redistribute
 *              it and/or modify it under the terms of the GNU General
 *              Public License as published by the Free Software
 *              Foundation; either version 2 of the License, or
 *              (at your option) any later version.
 *                This program is distributed in the hope that it will
 *              be useful, but WITHOUT ANY WARRANTY; without even the
 *              implied warranty of MERCHANTABILITY or FITNESS FOR A
 *              PARTICULAR PURPOSE. See the GNU General Public License
 *              for more details.
 *                You should have received a copy of the GNU General
 *              Public License along with this program; if not, write
 *              to the Free Software Foundation, Inc., 675 Mass Ave,
 *              Cambridge, MA 02139, USA.
 *
 *   AUTHORS: Favey Etienne   <favey@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            31-08-1997  --  created
 * $Log: tools.c,v $
 * Revision 1.7  2001/08/31 08:28:41  cocard
 * *** empty log message ***
 *
 * Revision 1.6  2001/08/23 08:35:46  cocard
 * *** empty log message ***
 *
 * Revision 1.5  2001/08/22 16:22:19  cocard
 * *** empty log message ***
 *
 * Revision 1.4  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.3  1998/10/30 14:22:07  favey
 * changed config.h to global.h
 *
 * Revision 1.2  1998/05/20 17:35:39  cocard
 * combMofN functions added (MC)
 *
 * Revision 1.1  1998/04/22 12:34:09  favey
 * Initial revision
 *
 * Revision 0.3  1997/12/24 12:12:50  favey
 * followup of V. 0.29
 *
 */

/*
 *  Standard Libraries
 */
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/*
 *  our special headers
 */

#include "global.h"
#include "errproc.h"
#include "tools.h"

/*----------------------------------------------------------------------*/
int arc_intersect (
  const double p1[3],
  const double p2[3],
  double na[2],
  double nb[2])
{
  /* Those variables are only needed locally */
  extern int errno;
  double p1x, p1y, s1, p2x, p2y, s2;
  double s12, t12, alpha, t1na, t1nb;
  double tmp1;

  p1x = p1[0];
  p1y = p1[1];
  s1  = p1[2];
  p2x = p2[0];
  p2y = p2[1];
  s2  = p2[2];

  /* calculation */
  s12 = sqrt ((p2x-p1x)*(p2x-p1x) + (p2y-p1y)*(p2y-p1y));

  /* check whether points are identical */
  if (s12 == 0.0) {
#ifdef WITH_ERRORS
    fprintf (stderr, "GGL Common Library --- tools\n");
    fprintf (stderr, "arc_intersect: points identical: %s\n",
	     strerror (ERANGE));
#endif
    errno = ERANGE;
    return EINVAL;
  }

  /* check whether arc cosine is greater or smaller than M_PI */
  if (p1y<=p2y)
    t12 = acos ((p2x-p1x)/s12);
  else 
    t12 = 2*M_PI - acos ((p2x-p1x)/s12);
  
  tmp1 = (s1*s1 + s12*s12 - s2*s2)/(2*s1*s12);
  /* check whether points are close enough that arc will intersect */
  if (fabs (tmp1) > 1) {
#ifdef WITH_ERRORS
    fprintf (stderr, "GGL Common Library --- tools\n");
    fprintf (stderr, "arc_intersect: points too far apart: %s\n",
	     strerror (EDOM));
#endif
    errno = EDOM;
    return EINVAL;
  }
  alpha = acos (tmp1);
      
  /* solution one */
  t1na = t12 + alpha;
  na[0] = p1x + s1 * cos(t1na);
  na[1] = p1y + s1 * sin(t1na);
  
  /* solution two */
  t1nb = t12 - alpha;
  nb[0] = p1x + s1 * cos(t1nb);
  nb[1] = p1y + s1 * sin(t1nb);

#ifdef LIB_DEBUG
  fprintf (stderr,
	   "arc_intersect: s12 | t12 | alpha | t1na | t1nb | na | nb\n");
  fprintf (stderr, "%g | %g | %g | %g | %g | x %g | y %g | x %g | y %g\n",
	   s12, t12/M_PI*180, alpha/M_PI*180, t1na/M_PI*180, t1nb/M_PI*180,
	   na[0], na[1], nb[0], nb[1]);
#endif
  return 0;
}

/*----------------------------------------------------------------------
 * routines for picking all subset of M integer from a set of N integer
 * example for M = 3 and N = 4
 * 1:   0  1  2
 * 2:   0  1  3
 * 3:   0  2  3
 * 4:   1  2  3
 *  routines
 *      checkMofN    
 *      combMofN_first
 *      combMofN_next
 *      combMofN_count
 *---------------------------------------------------------------------- */


/*----------------------------------------------------------------------
 * error routine used by all combMofN routines
 */
static void checkMofN(const int M, const int N, char * routine)
{
  if (M==0 || N==0 || M>N)
    { 
      fprintf (stderr, "GGL Common Library --- tools:\n");
      fprintf (stderr, "in: %s (M=%d/N=%d)\n", routine, M, N);
    }
}

/*----------------------------------------------------------------------
 * initializes the first combination. All following
 * are obtained by using combMofN_next
 * NO memory allocation for array comb (must be done outside)
 */
int combMofN_first(int *comb, const int M)
{
  int i;
#ifdef WITH_ERRORS
  checkMofN(M,M,"combMofN_init");
#endif
  for (i=0; i<M; i++) comb[i]=i;
  return 0;
}

/*----------------------------------------------------------------------
 * replaces comb by the next combination. returns 0 on success
 * and -1 after the last combination was reached. 
 * example :
 * combMofN_first(4,3) comb = 0 1 2  (return value =  0)
 * combMofN_next(4,3)  comb = 0 1 3  (return value =  0)
 * combMofN_next(4,3)  comb = 0 2 3  (return value =  0)
 * combMofN_next(4,3)  comb = 1 2 3  (return value =  0)
 * combMofN_next(4,3)  comb = 1 2 3  (return value = -1 -> no new combination returned)
 */
int combMofN_next (int *comb, const int M, const int N)
{
  int i=M;
  int j=N;

#ifdef WITH_ERRORS
  checkMofN(M,N,"combMofN_next");
#endif

  do {
    i--; 
    j--;
    if (i<0) return -1;     /* no valid new combination possible */
  } while (comb[i]==j);
  
  comb[i]++;
  for (j=i+1;j<M;j++) comb[j]=comb[j-1]+1; 
  return 0;
}

/*----------------------------------------------------------------------
 * returns the number of combinations (M of N) as a double 
 */
double combMofN_count(const int M, const int N)
{
  double count;
  int i, j, k;

#ifdef WITH_ERRORS
  checkMofN(M,N,"combMofN_count");
#endif
  
  if (M==1) return N;
  if (N==2) return 1;
  j =  N-M;
  j = M>j ? M : j;
  k = N - j;
  count=1;
  for (i=N; i>k; i--) count *= i;
  for (i=j; i>1; i--) count /= i;
  return  count;
}

/*----------------------------------------------------------------------
 * routines for determining all permutation within an array
 * example for count = 3
 * 1:   0  1  2
 * 2:   1  0  2
 * 3:   0  2  1
 * 4:   2  0  1
 * 5.   1  2  0
 * 6:   2  1  0
 *  routines
 *      permutation_init  
 *      permutation_next 
 *      permutation_free
 *      permutation_ctrlprint
 *
 * ( for good understanding of how to use this routines see
 *   the code of permutation_ctrlprint )
 *---------------------------------------------------------------------- */


/*----------------------------------------------------------------------
 * initialize memory for the array <sequence>,  and for an additional
 * control array <ctrlarray> which is needed in permutation_next
 * n is the dimension of the array
 * fill the array with the first sequence
 */
 int permutation_init(int n, int **sequence, int **ctrlarray)
{
  int i;

  if (((*sequence)  = (int *) malloc(sizeof(int)*n)) == NULL) return -1;
  if (((*ctrlarray) = (int *) malloc(sizeof(int)*(n+1))) == NULL) return -1;

  for (i=0; i<n; i++) (*sequence)[i]=i;
  memset((*ctrlarray), 0, sizeof(int)*n);
  (*ctrlarray)[n] = 1; /* some position index needed in permutation_next */
  return 0;
}

/*--------------------------------------------------------------------
 * permutation to the next sequence
 */
int permutation_next(int n, int * sequence, int * ctrlarray)
{
  int i, j, k;
  int pos = ctrlarray[n];

  while (pos < n)
    {
      i=sequence[pos];
      ctrlarray[i]=0;
      do i--; while (i >= 0 && ctrlarray[i]==1);
      if (i >= 0)
	{
	  sequence[pos]=i;
	  ctrlarray[i]=1;
	  j = k = 0;
	  while (j<pos)
	    {
	      if (ctrlarray[k] == 0)  sequence[j++]=k;  
	      k++;   
	    }
	  ctrlarray[n] = 1;
	  return 0;
	}
      pos++;
    }
  ctrlarray[n]=pos;
  return -1;
}
 
/*--------------------------------------------------------------------
 * returns the number of possible permutations
 */
long int permutation_count(int n)
{
  long int total = 1;
  int i;

  if (n < 1) return 0;
  for (i=2; i<=n; i++)
    total *= i;
  return total;
}

/*--------------------------------------------------------------------
 * free the memory allocated by permutation_init
 */
int permutation_free(int *sequence, int *ctrlarray)
{
  free(sequence);
  free(ctrlarray);
  return 0;
}

/*--------------------------------------------------------------------
 * routine for writing a complete set of all permutations 
 * of n integer in the range from 0..(n-1) to a file 
 * if addcounter == true the lines are numbered
 */
int permutation_ctrlprint(FILE *fp, int n, int addcounter)
{
  int *seq, *ctrl;
  long int total, count;
  char format[80];
  int i;

  /* some cosmetics for writing the numbering */
  count = total = permutation_count(n);  
  i=1;  while ((count = count/10) != 0) i++;
  sprintf(format,"%%0%dld/%%0%dld:",i, i);

  /* main loop */
  count=0;
  if (permutation_init(n, &seq, &ctrl) != 0) return -1;
  do
    {
      count++;
      if (addcounter) fprintf(fp,format,count,total);
      for (i=0; i<n; i++) fprintf(fp," %2d",seq[i]);
      fprintf(fp,"\n");
    }
  while (permutation_next(n, seq, ctrl) == 0);
  permutation_free(seq, ctrl);
  return 0;  
}
