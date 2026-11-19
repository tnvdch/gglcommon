/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/sed_mat3.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations on 3x3 matrix and 3x1 vectors
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: sed_mat3.c,v 1.3 2001/02/21 13:28:19 cocard Exp favey $
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
 *   AUTHORS: Cocard Marc     <cocard@geod.ethz.ch>
 *            Favey Etienne   <favey@geod.ethz.ch>
 *
 *   HISTORY: 
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/sed_mat3.c $
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.3  2001/02/21 13:28:19  cocard
 * *** empty log message ***
 *
 * Revision 1.2  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.1  1999/02/01 17:02:54  favey
 * Initial revision
 *
 * Revision 1.2  1998/10/30 14:22:07  favey
 * changed config.h to global.h
 *
 * Revision 1.1  1998/04/22 12:33:34  favey
 * Initial revision
 *
 * Revision 0.4  1998/01/19 15:14:03  favey
 * added function m3_write()
 *
 * Revision 0.3  1997/12/24 12:11:25  favey
 * followup of V. 0.29
 *
 *
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
#include "SED_Tmatcalc.h"
#include "SED_Tmat3.h"

int SED_Tm3_add ( SED_MATTYPE a[3][3],
		  SED_MATTYPE b[3][3],
		  SED_MATTYPE c[3][3])
{
  int i,j;
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      c[i][j] = a[i][j] + b[i][j];
  return 0;
}

int SED_Tm3_sub ( SED_MATTYPE a[3][3],
		  SED_MATTYPE b[3][3],
		  SED_MATTYPE c[3][3])
{
  int i,j;
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      c[i][j] = a[i][j] - b[i][j];
  return 0;
}

int SED_Tm3_unit (SED_MATTYPE a[3][3])
{
  int i,j;
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      a[i][j] = (i==j) ?   1.0 :   0.0;
  return 0;
}

int SED_Tm3_copy ( SED_MATTYPE a[3][3],
		   SED_MATTYPE c[3][3])
{
  int i,j;
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      c[i][j] = a[i][j];
  return 0;
}

int SED_Tm3_trans ( SED_MATTYPE a[3][3],
		    SED_MATTYPE c[3][3])
{
  SED_MATTYPE b[3][3];
  int i,j;
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      b[i][j] = a[j][i];
  SED_Tm3_copy (b, c);
 return 0;
}

int SED_Tm3_inv ( SED_MATTYPE a[3][3])
{
  struct SED_longTmatrix *m;
  int i,j,err;
  m = SED_Tm_init(3,3);
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      m->mat[i][j] = a[i][j];
  err = SED_Tm_inv (m);
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      a[i][j] = m->mat[i][j];
  SED_Tm_free (m);
  return err;
}

void SED_Tm3_write (
	FILE *fpout,
	 SED_MATTYPE a[3][3])
{
  int i,j;
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++)
      fprintf (fpout, "%12.8SED_bigTf ", a[i][j]);
    fprintf (fpout, "\n");
  }
  fprintf (fpout, "\n");
}

int SED_Tv3_add ( SED_MATTYPE a[3],
		  SED_MATTYPE b[3],
		  SED_MATTYPE c[3])
{
  int i;
  for (i=0;i<3;i++) c[i] = a[i] + b[i];
  return 0;
}

int SED_Tv3_sub ( SED_MATTYPE a[3],
		  SED_MATTYPE b[3],
		  SED_MATTYPE c[3])
{
  int i;
  for (i=0;i<3;i++) c[i] = a[i] - b[i];
  return 0;
}

int SED_Tv3_copy ( SED_MATTYPE a[3],
		   SED_MATTYPE c[3])
{
  int i;
  for (i=0;i<3;i++) c[i] = a[i];
  return 0;
}


int SED_Tv3_multx ( SED_MATTYPE a[3],
		    SED_MATTYPE x,
		    SED_MATTYPE c[3])
{
  int i;
  for (i=0;i<3;i++) c[i] = a[i] * x;
  return 0;
}

int SED_Tv3_norm ( SED_MATTYPE a[3],
		   SED_MATTYPE c[3])
{
  SED_MATTYPE n;
  n = (  1.0) / SED_Tv3_length(a);
  SED_Tv3_multx(a, n, c);
  return 0;
}

SED_MATTYPE SED_Tv3_length ( SED_MATTYPE a[3])
{
  return  (SED_MATTYPE) sqrt ((double) SED_Tv3_scalar(a,a));
}

SED_MATTYPE SED_Tv3_scalar ( SED_MATTYPE a[3],
			     SED_MATTYPE b[3])
{
  return   (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

int SED_Tv3_cross ( SED_MATTYPE a[3],
		    SED_MATTYPE b[3],
		    SED_MATTYPE c[3])
{
  SED_MATTYPE l[3];
  l[0] = a[1]*b[2] - a[2]*b[1];
  l[1] = a[2]*b[0] - a[0]*b[2];
  l[2] = a[0]*b[1] - a[1]*b[0];
  SED_Tv3_copy (l, c);
  return 0;
}

int SED_Tm3v3_mult ( SED_MATTYPE a[3][3],
		     SED_MATTYPE b[3],
		     SED_MATTYPE c[3])
{
  SED_MATTYPE l[3];
  int i,j;
  for (i=0;i<3;i++)
    {
      l[i] =   0.0;
      for (j=0;j<3;j++)
	l[i] += a[i][j] * b[j];
    }
  SED_Tv3_copy (l, c);
  return 0;
}

int SED_Tm3m3_mult ( SED_MATTYPE a[3][3],
		     SED_MATTYPE b[3][3],
		     SED_MATTYPE c[3][3])
{
  SED_MATTYPE l[3][3];
  int i,j,k;
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      {
	l[i][j] =   0.0;
	for (k=0;k<3;k++)
	l[i][j] += a[i][k] * b[k][j];
      }
  SED_Tm3_copy (l, c);
  return 0;
}

int SED_Trot_getmat ( int axis,
		      SED_MATTYPE angle,
		      SED_MATTYPE rotmat[3][3])
{
  SED_Tm3_unit (rotmat);
  switch (axis)
    {
    case 1:
      rotmat[1][1] = rotmat[2][2] = cos(angle);
      rotmat[2][1] = - (rotmat[1][2] = sin(angle));
      break;
    case 2:
      rotmat[0][0] = rotmat[2][2] = cos(angle);
      rotmat[0][2] = - (rotmat[2][0] = sin(angle));
      break;
    case 3:
      rotmat[0][0] = rotmat[1][1] = cos(angle);
      rotmat[1][0] = - (rotmat[0][1] = sin(angle));
      break;
    default: return -1;
    }
  return 0;
}

int SED_Trotate ( SED_MATTYPE x1[3],
		  int axis,
		  SED_MATTYPE angle,
		  SED_MATTYPE target[3])
{
  SED_MATTYPE x2[3];
   SED_MATTYPE co = cos(angle);
   SED_MATTYPE si = sin(angle);
  switch (axis)
    {
    case 1:
      x2[0] = x1[0];
      x2[1] =             x1[1]*co  +  x1[2]*si;
      x2[2] =            -x1[1]*si  +  x1[2]*co;
      break;
    case 2:
      x2[0] = x1[0]*co              -  x1[2]*si;
      x2[1] =             x1[1];
      x2[2] = x1[0]*si              +  x1[2]*co;
      break;
    case 3:
      x2[0] = x1[0]*co  +  x1[1]*si;
      x2[1] =-x1[0]*si  +  x1[1]*co;
      x2[2] =                          x1[2];
      break;
    default: return -1;
    }
  SED_Tv3_copy (x2, target);
  return 0;
}

int SED_Trotvec_to_mat ( SED_MATTYPE rotvec[3],
			 SED_MATTYPE angle,
			 SED_MATTYPE rotmat[3][3])
{
  SED_MATTYPE x = sin(angle);
  int i,j;
  rotmat[1][0] = -(rotmat[0][1] =  rotvec[2]*x);
  rotmat[0][2] = -(rotmat[2][0] =  rotvec[1]*x);
  rotmat[2][1] = -(rotmat[1][2] =  rotvec[0]*x);
  x = cos(angle);
  for (i=0;i<3;i++) rotmat[i][i] = x;
  x = 1 - x;
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      rotmat[i][j] += x*rotvec[i]*rotvec[j];
  return 0;
}

int SED_Trotmat_to_vec ( SED_MATTYPE rotmat[3][3],
			 SED_MATTYPE rotvec[3],
			 SED_MATTYPE *angle)
{
  int i;
  SED_MATTYPE x;

  *angle = 0;
  for (i=0;i<3;i++) *angle += rotmat[i][i];
  x = 2* sin(*angle = acos((*angle - 1.0)/2.0));
  if ( x  == 0.0 )
    return -1;
  rotvec[0] = (rotmat[1][2] - rotmat[2][1]) / x;
  rotvec[1] = (rotmat[2][0] - rotmat[0][2]) / x;
  rotvec[2] = (rotmat[0][1] - rotmat[1][0]) / x;
  return 0;
}
