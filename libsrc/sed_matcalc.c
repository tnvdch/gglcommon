/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/sed_matcalc.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations (parent file for i- l- and ordinary matcalc)
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: sed_matcalc.c,v 1.15 2001/08/31 08:28:41 cocard Exp favey $
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
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/sed_matcalc.c $
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.15  2001/08/31 08:28:41  cocard
 * *** empty log message ***
 *
 * Revision 1.14  2001/02/21 14:26:58  cocard
 * *** empty log message ***
 *
 * Revision 1.13  2000/05/12 09:22:53  favey
 * changed devel version to DEVEL_MAT_DEBUG
 *
 * Revision 1.12  2000/04/17 15:54:16  cocard
 * *** empty log message ***
 *
 * Revision 1.11  2000/02/25 12:57:58  ypeter
 * m_polynom_approx: modified function call, added possible weightedobs.
 * m_polynom_interpol: moved from winterpol to matcalc
 *
 * Revision 1.10  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.9  1999/08/18 08:58:38  cocard
 * Revision 1.8  1999/08/17 13:12:13  cocard
 * changes in functions min-max-absmax -- NULL pointer for row/col allowed
 * Revision 1.7  1999/07/29 13:22:28  cocard
 * $new function m_resize added (MC)
 * Revision 1.6  1999/06/03 10:09:13  cocard
 * Revision 1.5  1999/06/02 16:13:19  cocard
 * Revision 1.4  1999/04/19 09:51:06  cocard
 * new routine m_det added
 * Revision 1.3  1999/04/08 12:32:29  cocard
 * fclose added in m_readascii
 * Revision 1.2  1999/04/07 11:45:32  cocard
 * read/write binary debugged
 * Revision 1.1  1998/11/03 18:00:51  favey
 * Initial revision
 * Revision 1.8  1998/10/30 14:22:07  favey
 * changed config.h to global.h
 * Revision 1.7  1998/08/14 15:42:03  cocard
 * bugfix in m_mult (MC)
 * Revision 1.6  1998/07/01 14:30:54  cocard
 * m_separate_sym added
 * Revision 1.5  1998/06/24 09:47:32  cocard
 * conditional DEVEL_MAT_DEBUG removed for error routines m_error themselves
 * Revision 1.4  1998/06/23 16:06:10  favey
 * Revision 1.3  1998/06/08 11:59:18  favey
 * fixed bug in choleskyinv
 * Revision 1.2  1998/05/06 08:40:34  cocard
 * Revision 1.1  1998/04/22 12:33:18  favey
 * Revision 0.5  1998/01/07 15:43:06  favey
 * fixed bugs in m_multx and m_atpa1 (thanks to Mario Cerniar)
 * Revision 0.4  1998/01/05 15:54:35  favey
 * fixed allocation bug in m_copy
 * Revision 0.3  1997/12/24 12:11:33  favey
 * followup of V. 0.29
 *
 *   CAVEATS: contains error messages exiting main()
 *            when compiled with DEVEL_MAT_DEBUG defined
 *            returning EXIT_FAILURE
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
 *  our special header
 */
#include "global.h"
#include "errproc.h"
#include "SED_Tmatcalc.h"

/*****************************************************************/
void /*const*/ SED_Tm_WriteMatInfo (struct SED_longTmatrix *m, int i)
{
      fprintf (stderr, "SED_Tmatcalc: SED_longTmatrix %d <%s>\n", 
	       i, m->name==NULL?"NoName":m->name);
      fprintf (stderr, "SED_Tmatcalc: dimension : %d x %d\n",
	       m->row, m->col);
      fprintf (stderr, "SED_Tmatcalc: mem alloc : %d x %d\n",
	       m->dim_row, m->dim_col);  
}

/*****************************************************************/
void SED_Tm_error (
	const char *procname,
	struct SED_longTmatrix *m1, 
	struct SED_longTmatrix *m2,
	struct SED_longTmatrix *m3,
	int *row, 
	int *col)
{
  fprintf (stderr, "GGL Common Library -- SED_Tmatcalc\n");
  fprintf (stderr, "ERROR in subroutine: %s\n", procname);
  if (m1 != NULL) SED_Tm_WriteMatInfo (m1,1);
  if (m2 != NULL) SED_Tm_WriteMatInfo (m2,2);
  if (m3 != NULL) SED_Tm_WriteMatInfo (m3,3);
  if (row != NULL) fprintf (stderr, "SED_Tmatcalc: operation on row %d\n", *row);
  if (col != NULL) fprintf (stderr, "SED_Tmatcalc: operation on col %d\n", *col);
  fprintf (stderr,"\n");
}



/*****************************************************************/
struct SED_longTmatrix * /*const*/ SED_Tm_get (int row, int col)
{
  struct SED_longTmatrix *m;
  SED_MATTYPE *ptmp;
  int i, j;

  m = (struct SED_longTmatrix *) malloc (sizeof(struct SED_longTmatrix));
  if ( m == NULL )
  {
 #ifdef DEVEL_MAT_DEBUG
       SED_Tm_error ("SED_Tm_init: virtual memory exhausted",
 	       NULL, NULL, NULL, NULL, NULL);
 #endif
       return NULL;
  }

  m->fast = (SED_MATTYPE *) malloc (row*col*sizeof(SED_MATTYPE));
  if ( m->fast == NULL )
  {
  #ifdef DEVEL_MAT_DEBUG
	SED_Tm_error("SED_Tm_init: virtual memory exhausted",
  	   NULL, NULL, NULL, NULL, NULL);
  #endif
	free(m);
  	return NULL;
  }
  m->mat = (void *) malloc (row*sizeof(SED_MATTYPE *));
  if ( m->mat == NULL )
  {
  #ifdef DEVEL_MAT_DEBUG
      SED_Tm_error ("SED_Tm_init: virtual memory exhausted",
  	       NULL, NULL, NULL, NULL, NULL);
  #endif
	free(m->fast);
	free(m);
	return NULL;
  }
  m->row = m->dim_row = row;
  m->col = m->dim_col = col;
  m->name = NULL;

  ptmp = m->fast;
  for (i=0;i<row;i++)
    {
      m->mat[i] = ptmp;
      for (j=0;j<col;j++) ptmp++;
    }

  return m;
}

/*****************************************************************/
struct SED_longTmatrix * SED_Tm_init (int row, int col)
{
  struct SED_longTmatrix *m;
  m = SED_Tm_get (row, col);
  if (m != NULL) SED_Tm_zero(m);
  return m;
}

/*****************************************************************/
int SED_Tm_name (struct SED_longTmatrix *m, const char *newname)
{
  if (m->name!=NULL) free(m->name);
  if (strlen(newname)>0) m->name = strdup(newname);
  return 0;
}

/*****************************************************************/
struct SED_longTmatrix * SED_Tm_init_name(int row, int col, char *newname)
{
  struct SED_longTmatrix *m;
  if ((m = SED_Tm_init (row, col)) == NULL) return NULL;
  SED_Tm_name (m, newname);
  return m;
}

/*****************************************************************/
void SED_Tm_free (struct SED_longTmatrix *m)
{
  if (m->name!=NULL) free (m->name);
  free (m->mat);
  free (m->fast);
  free (m);
}

/*****************************************************************/
int SED_Tm_resize (struct SED_longTmatrix *m, 
		   int row, int col, int set_zero)
{
#ifdef DEVEL_MAT_DEBUG
  if (row > m->dim_row || col > m->dim_col) 
    {
      SED_Tm_error ("SED_Tm_resize", m, NULL, NULL, &row, &col);
      return -1; 
    }
#endif
  m->row=row;
  m->col=col;
  if (set_zero) SED_Tm_zero(m);
  return 0;
}

/*****************************************************************/
int SED_Tm_fill (struct SED_longTmatrix *m,
		 int row, int col, SED_MATTYPE value)
{
#ifdef DEVEL_MAT_DEBUG
  if ( (row>m->row) || (col>m->col)
       || (row<0) || (col<0) )
    {
      SED_Tm_error ("SED_Tm_fill", m, NULL, NULL, &row, &col);
      return -1;
    }
#endif

  m->mat[row][col] = value;
  return 0;
}

/*****************************************************************/
int SED_Tm_fill_row (struct SED_longTmatrix *m,
		     int row, SED_MATTYPE value)
{
  int i;

#ifdef DEVEL_MAT_DEBUG
  if ( (row>m->row)
       || (row<0) )
    {
      SED_Tm_error ("SED_Tm_fill_row", m, NULL, NULL, &row, NULL);
      return -1;
    }
#endif

  for (i=0;i<m->col;i++) m->mat[row][i] = value;
  return 0;
}

/*****************************************************************/
int SED_Tm_fill_col (struct SED_longTmatrix *m,
		     int col, SED_MATTYPE value)
{
  int i;

#ifdef DEVEL_MAT_DEBUG
  if ( (col>m->col)
       || (col<0) )
    {
      SED_Tm_error ("SED_Tm_fill_col", m, NULL, NULL, NULL, &col);
      return -1;
    }
#endif

  for (i=0;i<m->row;i++) m->mat[i][col] = value;
  return 0;
}

/*****************************************************************/
int SED_Tm_fillsym (struct SED_longTmatrix *m)
{
  int i, j;

#ifdef DEVEL_MAT_DEBUG
  if (m->row != m->col) 
    {
      SED_Tm_error ("SED_Tm_fillsym", m, NULL, NULL, NULL, NULL);
      return -1;
    }
#endif
  for (i=0;i<m->row;i++)
    for (j=0;j<i;j++)
      m->mat[i][j] = m->mat[j][i];
  return 0;
}

/*****************************************************************/
int SED_Tm_zero (struct SED_longTmatrix *m)
{
int i;

#ifdef DEVEL_MAT_DEBUG
  if (m->fast == NULL || m->mat == NULL)
    SED_Tm_error ("SED_Tm_zero\nplease initalize your matrix first",
	     m, NULL, NULL, NULL, NULL);
#endif
  if (m->dim_row == m->row && m->dim_col == m->col)
    memset (m->fast, 0, m->dim_row*m->dim_col*sizeof(SED_MATTYPE));
  else
    for (i=0; i<m->row; i++)
      memset (m->mat[i], 0, m->col*sizeof(SED_MATTYPE));  
  return 0;
}

/*****************************************************************/
int SED_Tm_zero_row (struct SED_longTmatrix *m,
	        int row)
{
#ifdef DEVEL_MAT_DEBUG
  if (m->dim_row<row)
    {
      SED_Tm_error ("SED_Tm_zero_row", m, NULL, NULL, &row, NULL);
      return -1;
    }
#endif

  memset (m->mat[row], 0, m->col*sizeof(SED_MATTYPE));
  return 0;
}

/*****************************************************************/
int SED_Tm_unit (struct SED_longTmatrix *m)
{
  int i;
#ifdef DEVEL_MAT_DEBUG
  if (m->row != m->col) 
    {
      SED_Tm_error ("SED_Tm_unit\nmatrix not diagonal", m, NULL, NULL, NULL, NULL);
      return -1;
    }
#endif

  SED_Tm_zero(m);
  for (i=0;i<m->col;i++) m->mat[i][i] = (SED_MATTYPE) 1.0 ;
  return 0;
}

/*****************************************************************/
int SED_Tm_set (struct SED_longTmatrix *m,
	   SED_MATTYPE value)
{
  int i, j;

  for (i=0;i<m->row;i++)
    for (j=0;j<m->col;j++)
      m->mat[i][j] = value;
  return 0;
}

/*****************************************************************/
int SED_Tm_random (struct SED_longTmatrix *m, SED_MATTYPE fac)
{
  int i, j;

  for (i=0;i<m->row;i++)
    for (j=0;j<m->col;j++)
      m->mat[i][j] = (SED_MATTYPE) fac*((double) rand()/ RAND_MAX);  
  return 0;
}

/*****************************************************************/
int SED_Tm_copy (struct SED_longTmatrix *source,
		 struct SED_longTmatrix *target)
{
  int i, j;

#ifdef DEVEL_MAT_DEBUG
  if ( (target->dim_row < source->row) || (target->dim_col < source->col) )
    {
      SED_Tm_error ("SED_Tm_copy", source, target, NULL, NULL, NULL);
      return -1;
    }
#endif

  target->row = source->row;
  target->col = source->col;
  if (target->name != NULL) free (target->name);
  target->name = (source->name == NULL) ? NULL : strdup (source->name);
  for (i=0;i<source->row;i++)
    for (j=0;j<source->col;j++)
      target->mat[i][j] = source->mat[i][j];  
  return 0;
}

/*****************************************************************/
int SED_Tm_pushup (struct SED_longTmatrix *m)
{
  if (m->row < m->dim_row)
    {
      memmove (m->mat[m->row], m->mat[0], m->dim_col*sizeof(SED_MATTYPE));
      memmove (m->mat[0], m->mat[1], m->row*m->dim_col*sizeof(SED_MATTYPE));
    }
  else
    {
      SED_MATTYPE *ptmp;

      ptmp = (SED_MATTYPE *) malloc (m->dim_col*sizeof(SED_MATTYPE));
      memcpy (ptmp, m->mat[0], m->dim_col*sizeof(SED_MATTYPE));
      memmove (m->mat[0], m->mat[1], (m->row-1)*m->dim_col*sizeof(SED_MATTYPE));
      memcpy (m->mat[m->row-1], ptmp, m->dim_col*sizeof(SED_MATTYPE));
      free (ptmp);
    }
  return 0;
}

/*****************************************************************/
int SED_Tm_pushdown (struct SED_longTmatrix *m)
{
  if (m->row < m->dim_row)
    {
      memmove (m->mat[1], m->mat[0], m->row*m->dim_col*sizeof(SED_MATTYPE));
      memmove (m->mat[0], m->mat[m->row], m->dim_col*sizeof(SED_MATTYPE));
    }
  else
    {
      SED_MATTYPE *ptmp;

      ptmp = (SED_MATTYPE *) malloc (m->dim_col*sizeof(SED_MATTYPE));
      memcpy (ptmp, m->mat[m->row-1], m->dim_col*sizeof(SED_MATTYPE));
      memmove (m->mat[1], m->mat[0], (m->row-1)*m->dim_col*sizeof(SED_MATTYPE));
      memcpy (m->mat[0], ptmp, m->dim_col*sizeof(SED_MATTYPE));
      free (ptmp);
    }
  return 0;
}

/*****************************************************************/
int SED_Tm_max_col (struct SED_longTmatrix *m,
	       const int col,
	       SED_MATTYPE *max,
	       int *row)
{
  int i, r;

  *max = m->mat[r=0][col];
  for (i=1;i<m->row;i++) 
    *max = ((m->mat[i][col]>*max) ? m->mat[r=i][col] : *max);
  if (row != NULL) (*row)=r;
  return 0;
}

/*****************************************************************/

int SED_Tm_min_col (struct SED_longTmatrix *m,
	       const int col,
	       SED_MATTYPE *min,
	       int *row)
{
  int i, r;

  *min = m->mat[r=0][col];
  for (i=1;i<m->row;i++) 
    *min = ((m->mat[i][col]<*min) ? m->mat[r=i][col] : *min);
  if (row != NULL) (*row)=r;
  return 0;
}

/*****************************************************************
 *  the following funcs could be faster using m->fast, if
 *  logical and physical size are equal
 */

int SED_Tm_max (struct SED_longTmatrix *m,
	   SED_MATTYPE *max,
	   int *row,
	   int *col)
{
  int i, j, r, c;

  *max = m->mat[r=0][c=0];
  
  for (i=0;i<m->row;i++)
    for (j=0;j<m->col;j++)
      *max = ((m->mat[i][j]>*max) ? m->mat[r=i][c=j] : *max);
  if (row != NULL) (*row)=r;
  if (col != NULL) (*col)=c;  
  return 0;
}

/*****************************************************************/

int SED_Tm_min (struct SED_longTmatrix *m,
	   SED_MATTYPE *min,
	   int *row,
	   int *col)
{
  int i, j;
  int r, c;

  *min = m->mat[r=0][c=0];
  
  for (i=0;i<m->row;i++)
    for (j=0;j<m->col;j++)
      *min = ((m->mat[i][j]<*min) ? m->mat[r=i][c=j] : *min);
  if (row != NULL) (*row)=r;
  if (col != NULL) (*col)=c;  
  return 0;
}

/*****************************************************************/

int SED_Tm_absmax (struct SED_longTmatrix *m,
	   SED_MATTYPE *max,
	   int *row,
	   int *col)
{
  int i, j;
  int r, c;

  *max = m->mat[r=0][c=0];
  
  for (i=0;i<m->row;i++)
    for (j=0;j<m->col;j++)
      *max = (((SED_MATTYPE) fabs((double) m->mat[i][j])>*max) 
	      ? (SED_MATTYPE) fabs((double) m->mat[r=i][c=j]) : *max);
  if (row != NULL) (*row)=r;
  if (col != NULL) (*col)=c;  
  return 0;
}

/*****************************************************************/

void SED_Tm_writename (FILE *fpout,
		  struct SED_longTmatrix *m,
		  int tr)
{
   if (m->name != NULL) 
   {
     fprintf(fpout,"matrix: %s", m->name);
     if (tr != FALSE) fprintf(fpout," (transposed)");
       fprintf(fpout,"\n");
   }
}

/*****************************************************************/

int SED_Tm_write (FILE *fpout,
	     struct SED_longTmatrix *m,
	     const int fmt1,
	     const int fmt2,
	     int n,
             const int tr,
	     const int title)
{
  int i, j;
  int row, col;
  char mask[256], titlemask[256];

  row = ((tr != FALSE) ? m->col : m->row);
  col = ((tr != FALSE) ? m->row : m->col);
  if (n==0 || n>col) n=col;
 
  if (title != FALSE)
    {
      SED_Tm_writename(fpout, m, tr);
      sprintf (titlemask, " %%%dd", fmt1);
      for (j=0; j<col; j++)
	{
          if (j%n==0) fprintf (fpout, "   ");
	  fprintf (fpout, titlemask, (j+1));
	  if ((j+1) % n == 0 || (j+1) == col) fprintf (fpout, "\n");
	}
    }

  sprintf (mask, " %%%d.%df", fmt1, fmt2);
  for (i=0; i<row; i++)
    {
      for (j=0; j<col; j++)
        {
          if (title != FALSE && j%n==0) fprintf (fpout,"%3d", (i+1));
          if (tr != FALSE) fprintf (fpout, mask, m->mat[j][i]);
           else fprintf (fpout, mask, m->mat[i][j]);
	  if ((j+1) % n  == 0 || (j+1)==col) fprintf(fpout, "\n");
        }
    }
  return 0;
}
/*****************************************************************/

int SED_Tm_writefmt (FILE *fpout,
		struct SED_longTmatrix *m,
		const char * fmt,
		int n,
		const int tr,
		const int title)
{
  int i, j;
  int row, col;
  char mask[256];

  row = (tr != FALSE) ? m->col : m->row;
  col = (tr != FALSE) ? m->row : m->col;
  if (n==0 || n>col) n=col;
 
  if (title != FALSE)
    {
      sprintf (mask, fmt, m->mat[0][0]); 
      sprintf (mask, "%%%dd", (int) strlen(mask));
      SED_Tm_writename(fpout, m, tr); 
      for (j=0; j<col; j++)
	{
          if (j%n==0) fprintf (fpout, "   ");
	  fprintf (fpout, mask, (j+1));
	  if ((j+1) % n == 0 || (j+1) == col) fprintf (fpout, "\n");
	}
    }

  for (i=0; i<row; i++)
    {
      for (j=0; j<col; j++)
        {
          if (title != FALSE && j%n==0) fprintf (fpout,"%3d", (i+1));
          if (tr != FALSE) fprintf (fpout, fmt, m->mat[j][i]);
           else fprintf (fpout, fmt, m->mat[i][j]);
	  if ((j+1) % n  == 0 || (j+1)==col) fprintf(fpout, "\n");
        }
    }
  return 0;
}

/*****************************************************************/

int SED_Tm_write_filled (FILE *fpout,
	            struct SED_longTmatrix *m,
                    const int tr,
	            SED_MATTYPE nearlyzero)
{
  int i, j;
  int imax, jmax;

  imax = ((tr != FALSE) ? m->col : m->row);
  jmax = ((tr != FALSE) ? m->row : m->col);

  SED_Tm_writename(fpout, m, tr);
  if (imax>10)
    {
      fprintf(fpout,"    ");
      for (j=0; j<jmax;j+=10) fprintf(fpout,"%d         ",j/10);
     }
  fprintf(fpout,"\n    ");
  for (j=0;j<jmax;j++) fprintf(fpout,"%d", j % 10);
  fprintf(fpout,"\n");

  for (i=0;i < imax;i++)
    {
      fprintf(fpout,"%3d ",i);
      for (j=0;j < jmax;j++)
        {
          if (fabs(m->mat[i][j]) < nearlyzero)
            fprintf (fpout, "0");
          else
            fprintf (fpout, "*");
        }
      fprintf(fpout, "\n");
    }
  return 0;
}

/*****************************************************************/

int SED_Tm_writebin (FILE *fpout,
		struct SED_longTmatrix *m)
{
  int i;

  fprintf (fpout, "MATRIX %d %d %d %s\n", m->row, m->col, 
	   (int) sizeof(SED_MATTYPE), m->name==NULL?"NONAME":m->name);
  for (i=0;i<m->row;i++)
    fwrite (m->mat[i], sizeof(SED_MATTYPE), m->col, fpout);
  fprintf(fpout, "\n");
  return 0;
}

/*****************************************************************/
int SED_Tm_writebin2 (char * fname, struct SED_longTmatrix *m)
{
  FILE *fp;
  int returnvalue;

  if ( !(fp = fopen (fname, "w")) ) 
    {
#ifdef DEVEL_MAT_DEBUG
      char errmsg[256];
      sprintf (errmsg, " Could not open file `%s': %s\n",
	       fname,
	       strerror (errno) );
      SED_Tm_error (errmsg, NULL, NULL, NULL, NULL, NULL);
#endif
      return -1;
    }

  returnvalue=SED_Tm_writebin(fp, m);
  fclose(fp);
  return returnvalue;
}
/*****************************************************************/

int SED_Tm_readbin (FILE *fpin,
	       struct SED_longTmatrix **m)
{
  char ident[256];
  char name[256];
  int i;
  int size, row, col;

  memset(name,  0, 256);
  memset(ident, 0, 256);

  if (fscanf (fpin, "%s %d %d %d %s\n", ident, &row, &col, &size, name) != 5) return -1;
  if (strcmp(ident,"MATRIX") != 0) return -2;
  if (size != sizeof(SED_MATTYPE)) return -3;

  if (*m == NULL) /* initialize matrix */
    {
      if ((*m = SED_Tm_get (row, col)) == NULL) return -4;
    }
  else /* check dimensions */
    {
      if ((row<=(*m)->dim_row) && (col<=(*m)->dim_col))
	{
	  (*m)->row = row;
	  (*m)->col = col;
	}
      else
	return -5;
    }
   
  /* fill in name of matrix if present */
  if ((*m)->name != NULL) free((*m)->name);
  if (strcmp(name,"NONAME")==0) (*m)->name=NULL;  else (*m)->name=strdup(name);
    
  /* get data from file */
  for (i=0; i<(*m)->row; i++)
    if (fread((*m)->mat[i], sizeof(SED_MATTYPE), (*m)->col, fpin) != (*m)->col) return -6;
  
  return 0;
}

/*****************************************************************/
int SED_Tm_readbin2 (char * fname, struct SED_longTmatrix **m)
{
  FILE *fp;
  int returnvalue;

  if ( !(fp = fopen (fname, "r")) ) 
    {
#ifdef DEVEL_MAT_DEBUG
      char errmsg[256];
      sprintf (errmsg, " Could not open file `%s': %s\n",
	       fname,
	       strerror (errno) );
      SED_Tm_error (errmsg, NULL, NULL, NULL, NULL, NULL);
#endif
      return -1;
    }

  returnvalue=SED_Tm_readbin(fp, m);
  fclose(fp);
  return returnvalue;
}

/*****************************************************************/
int SED_Tm_readascii( char *fname, 
		 struct SED_longTmatrix **m, 
		 int row2read,
		 int col2read)
{
  FILE *fpin;
  double value;
  long int valcount=0; 
  int i,j;

  if ( !(fpin = fopen (fname, "r")) ) 
    {
#ifdef DEVEL_MAT_DEBUG
      char errmsg[256];
      sprintf (errmsg, " Could not open file `%s': %s\n",
	       fname,
	       strerror (errno) );
      SED_Tm_error (errmsg, NULL, NULL, NULL, NULL, NULL);
#endif
      return -1;
    }

  if (row2read ==0 && col2read==0)
    {
#ifdef DEVEL_MAT_DEBUG
       SED_Tm_error("error in SED_Tm_readascii\n (row2read = 0 and col2read = 0\n",
	NULL, NULL, NULL, NULL, NULL);       
#endif
      fclose(fpin);
      return -1;
    }
  
  if (row2read ==0 || col2read == 0) 
    {
      while (fscanf (fpin, " %lg", &value) != EOF)  /* count  values on file */
	{
	  valcount++;  
	}
      rewind(fpin); 

      if (valcount ==0)
      {
    	  fclose(fpin);
    	  return -1;
      }

      if (row2read != 0) 
	col2read = valcount/row2read; 	
      else 
	row2read = valcount/col2read; 
    }


  if ((*m)==NULL)   /* initialize matrix */
  {
		if (((*m) = SED_Tm_init (row2read, col2read)) == NULL)
		{
			fclose(fpin);
			return -1;
		}
  }
  else  /* check dimensions */
    {
      if (col2read > (*m)->dim_col && row2read > (*m)->dim_row)
      {
    	  fclose(fpin);
    	  return -1;
      }
      (*m)->col = col2read;
    }

  i = j = 0;
  while (fscanf (fpin, " %lg", &value) != EOF && i<(*m)->dim_row)
    {
      (*m)->mat[i][j++]=(SED_MATTYPE) value;
      if (j>=col2read) { i++; j=0; }
    }
  (*m)->row = i;
  fclose(fpin);
  return 0;
}

/*****************************************************************/

int SED_Tm_add (struct SED_longTmatrix *m1,
	   struct SED_longTmatrix *m2,
	   struct SED_longTmatrix *m_res)
{
  int i, j;

#ifdef DEVEL_MAT_DEBUG
  if ( (m1->row != m2->row)
       || (m1->col != m2->col)
       || (m_res->dim_row < m1->row)
       || (m_res->dim_col < m1->col) )
    {
      SED_Tm_error ("SED_Tm_add", m1, m2, m_res, NULL, NULL);
      return -1;
    }
#endif
  m_res->row = m1->row;
  m_res->col = m1->col;
  for (i=0;i<m1->row;i++)
    for (j=0;j<m1->col;j++)
      m_res->mat[i][j] = m1->mat[i][j] + m2->mat[i][j];
  return 0;
}

/*****************************************************************/

int SED_Tm_sub (struct SED_longTmatrix *m1,
	   struct SED_longTmatrix *m2,
	   struct SED_longTmatrix *m_res)
{
  int i, j;

#ifdef DEVEL_MAT_DEBUG
  if ( (m1->row != m2->row)
       || (m1->col != m2->col)
       || (m_res->dim_row < m1->row)
       || (m_res->dim_col < m1->col) )
    {
      SED_Tm_error ("SED_Tm_sub", m1, m2, m_res, NULL, NULL);
      return -1;
    }
#endif
  m_res->row = m1->row;
  m_res->col = m1->col;
  for (i=0;i<m1->row;i++)
    for (j=0;j<m1->col;j++)
      m_res->mat[i][j] = m1->mat[i][j] - m2->mat[i][j];
  return 0;
}

/*****************************************************************/

static void SED_Tmult____ (struct SED_longTmatrix *h, 
	      struct SED_longTmatrix *ma, 
	      struct SED_longTmatrix *mb,
	      const int m1row,
	      const int m2col,
	      const int m1col)
{
  int i,j,k;
  for  (i=0;i<m1row;i++)
    for (j=0;j<m2col;j++)
      {
	h->mat[i][j] = (SED_MATTYPE) 0.0;
	for (k=0;k<m1col;k++)
	  h->mat[i][j] += (ma->mat[i][k])*(mb->mat[k][j]);
      }
}

static void SED_Tmult_T__ (struct SED_longTmatrix *h, 
	      struct SED_longTmatrix *ma, 
	      struct SED_longTmatrix *mb,
	      const int m1row,
	      const int m2col,
	      const int m1col)
{
  int i,j,k;
  for  (i=0;i<m1row;i++)
    for (j=0;j<m2col;j++)
      {
	h->mat[i][j] = (SED_MATTYPE) 0.0;
	for (k=0;k<m1col;k++)
	  h->mat[i][j] += (ma->mat[k][i])*(mb->mat[k][j]);
      }
}

static void SED_Tmult___T (struct SED_longTmatrix *h, 
	      struct SED_longTmatrix *ma, 
	      struct SED_longTmatrix *mb,
	      const int m1row,
	      const int m2col,
	      const int m1col)
{
  int i,j,k;
  for  (i=0;i<m1row;i++)
    for (j=0;j<m2col;j++)
      {
	h->mat[i][j] = (SED_MATTYPE) 0.0;
	for (k=0;k<m1col;k++)
	  h->mat[i][j] += (ma->mat[i][k])*(mb->mat[j][k]);
      }
}

static void SED_Tmult_T_T (struct SED_longTmatrix *h, 
	      struct SED_longTmatrix *ma, 
	      struct SED_longTmatrix *mb,
	      const int m1row,
	      const int m2col,
	      const int m1col)
{
  int i,j,k;
  for  (i=0;i<m1row;i++)
    for (j=0;j<m2col;j++)
      {
	h->mat[i][j] = (SED_MATTYPE) 0.0;
	for (k=0;k<m1col;k++)
	  h->mat[i][j] += (ma->mat[k][i])*(mb->mat[j][k]);
      }
}


int SED_Tm_mult (struct SED_longTmatrix *m1,
	    const int m1tr,
	    struct SED_longTmatrix *m2,
	    const int m2tr,
	    struct SED_longTmatrix *m_res)
{
  int m1row, m1col, m2row, m2col;
  struct SED_longTmatrix *m_tmp;

  if (m1tr != FALSE)
    {
      m1row = m1->col; m1col = m1->row;
    }
  else
    {
      m1row = m1->row; m1col = m1->col;
    }
  if (m2tr != FALSE)
    {
      m2row = m2->col; m2col = m2->row;
    }
  else
    {
      m2row = m2->row; m2col = m2->col;
    }
#ifdef DEVEL_MAT_DEBUG
  if ( (m1col!=m2row) 
       || (m_res->dim_row<m1row) 
       || (m_res->dim_col<m2col) )
    {
      SED_Tm_error("SED_Tm_mult", m1, m2, m_res, NULL, NULL);
    }
#endif
  
  if (m_res == m1 || m_res == m2) m_tmp = SED_Tm_get (m1row, m2col);
  else m_tmp = m_res;
  if (m_tmp == NULL) return -1;

  if (m1tr != FALSE)
    {
      if (m2tr != FALSE)
	SED_Tmult_T_T (m_tmp, m1, m2, m1row, m2col, m1col);
      else
	SED_Tmult_T__ (m_tmp, m1, m2, m1row, m2col, m1col);
    }
  else
    {
      if (m2tr != FALSE)
	SED_Tmult___T (m_tmp, m1, m2, m1row, m2col, m1col);
      else
	SED_Tmult____ (m_tmp, m1, m2, m1row, m2col, m1col);
    }

  if (m_res == m1 || m_res == m2)
    {
      SED_Tm_copy (m_tmp, m_res);
      SED_Tm_free (m_tmp);
    }
  else
    {
      m_res->row=m1row;
      m_res->col=m2col;
    }
  
  return 0;
}

/*****************************************************************/

int SED_Tm_multx (struct SED_longTmatrix *m,
	     SED_MATTYPE scalar,
	     struct SED_longTmatrix *m_res)
{
  int i, j;

#ifdef DEVEL_MAT_DEBUG
  if ( (m_res->dim_row<m->row) 
       || (m_res->dim_col<m->col) )
    SED_Tm_error("SED_Tm_multx", m, m_res, NULL, NULL, NULL);
#endif

  if (m_res == NULL) m_res = SED_Tm_get (m->row, m->col);
  if (m_res == NULL) return -1;

  m_res->col = m->col;
  m_res->row = m->row;
  for (i=0;i<m->row;i++)
    for (j=0;j<m->col;j++)
      m_res->mat[i][j] = m->mat[i][j] * scalar;
  return 0;
}

/*****************************************************************/

SED_MATTYPE SED_Tm_trace (struct SED_longTmatrix *m)
{
  int i;
  SED_MATTYPE sp;

#ifdef DEVEL_MAT_DEBUG
  if (m->col!=m->row)
    SED_Tm_error("SED_Tm_trace", m, NULL, NULL, NULL, NULL);
#endif

  sp = (SED_MATTYPE) 0.0;
  for (i=0;i<m->row;i++) sp += m->mat[i][i];
  return sp;
}

/*****************************************************************/

int SED_Tm_trans (struct SED_longTmatrix *m,
	     struct SED_longTmatrix *m_res)
{
  int i, j;
  struct SED_longTmatrix *h;

#ifdef DEVEL_MAT_DEBUG
  if ( (m_res->dim_row<m->col) || (m_res->dim_col<m->row) )
    SED_Tm_error("SED_Tm_trans", m, NULL, NULL, NULL, NULL);
#endif

  h = m_res;
  if (m_res == m) h = SED_Tm_get (m_res->dim_row, m_res->dim_col);
  h->col = m->row;
  h->row = m->col;

  for (i=0;i<m->row;i++)
    for (j=0;j<m->col;j++)
      h->mat[j][i] = m->mat[i][j];
  if (m_res == m)
    {
      SED_Tm_copy(h, m_res);
      SED_Tm_free(h);
    }
  return 0;
}

/*****************************************************************/

int SED_Tm_atpa1 (struct SED_longTmatrix *a,
	     struct SED_longTmatrix *p,
	     struct SED_longTmatrix *atpa)
{
  int i,j,k;

#ifdef DEVEL_MAT_DEBUG
  if ( (a->row!=p->col)
       || (atpa->dim_col<a->col)
       || (atpa->dim_row<a->col) )
    SED_Tm_error("SED_Tm_atpa1", a, p, atpa, NULL, NULL);
#endif

  atpa->row = a->col;
  atpa->col = a->col;
  for (i=0;i<atpa->row;i++)
    for (j=0;j<atpa->col;j++)
      {
	atpa->mat[i][j] = (SED_MATTYPE) 0.0;
	for (k=0;k<a->row;k++)
	  atpa->mat[i][j] += (a->mat[k][i])*(a->mat[k][j])*(p->mat[0][k]);
      }
  SED_Tm_fillsym (atpa);
  return 0;
}

/*****************************************************************/

int SED_Tm_det (struct SED_longTmatrix *m, SED_MATTYPE *det)
{
  struct   SED_longTmatrix *h;
  int dim, line;
  int i, j, k;
  SED_MATTYPE x;
 
  if (m->row != m->col)
    {
#ifdef DEVEL_MAT_DEBUG
      SED_Tm_error ("SED_Tm_det", m, NULL, NULL, &(m->row), &(m->col));
#endif
      return -1;
    }

  dim=m->row;
  if ((h = SED_Tm_init(dim, dim)) == NULL) return -1; 
  SED_Tm_copy(m, h);  

  *det = 1;  
  for (i=0; i<dim; i++)  
    {
      if (((SED_MATTYPE) fabs ((double) h->mat[i][i])) < NEARLYZERO)
	{
	  line=i;
	  while (line<dim && ((SED_MATTYPE) fabs((double)  h->mat[line][i])) < NEARLYZERO) line++;	 
	  
	  if (line==dim) 
	    {
	      *det=0;
	      SED_Tm_free(h);
	      return 0;
	    }
	  for (j=i; j< dim; j++) 
	    {
	      x=h->mat[i][j]; 
	      h->mat[i][j]= h->mat[line][j]; 
	      h->mat[line][j]=x;
	    }
	}
      
      *det *= h->mat[i][i];
      for (j=i+1; j<dim; j++)
	{
	  x=h->mat[j][i]/h->mat[i][i];
	  for (k=i; k<dim; k++) h->mat[j][k]-=x*h->mat[i][k];	      
	}
    }
  SED_Tm_free(h);
  return 0;
}

/*****************************************************************/

int SED_Tm_pivot (struct SED_longTmatrix *m,
	     int row,
	     int col)
{
  SED_MATTYPE pivot, ai;
  int i, j;

#ifdef DEVEL_MAT_DEBUG
  if ( (row>=m->row) || (col>=m->col) )
    SED_Tm_error ("SED_Tm_pivot", m, NULL, NULL, &row, &col);
#endif

  pivot = m->mat[row][col];
  for (i=0;i<m->row;i++)
    if (i!=row)
      {
	ai = m->mat[i][col] / pivot;
	for (j=0;j<m->col;j++)
	  m->mat[i][j] = m->mat[i][j] - ai*(m->mat[row][j]);
	m->mat[i][col] = ai;
      }
  for (j=0;j<m->col;j++)
    m->mat[row][j] = -m->mat[row][j] / pivot;
  m->mat[row][col] = ((SED_MATTYPE) 1.0) / pivot;
  return 0;
}

/*****************************************************************/

int SED_Tm_inv (struct SED_longTmatrix *m)
{
  int nr=0;
  int ef;
  SED_MATTYPE pivot, ai, max, x;
  int *pivotcontrol;
  int stop;
  SED_MATTYPE nzero;
  int i, j;

#ifdef DEVEL_MAT_DEBUG
  if (m->row!=m->col)
    SED_Tm_error("SED_Tm_inv", m, NULL, NULL, NULL, NULL);
#endif

  pivotcontrol = (int *) malloc (m->col*sizeof(int));

  nzero = SED_Tm_trace (m);
  if (nzero!= (SED_MATTYPE) 0.0)
    nzero *= NEARLYZERO*m->row;
  stop = FALSE;
  for (i=0;i<m->col;i++) pivotcontrol[i] = TRUE;
  ef = 0;
  do
    {
      max = (SED_MATTYPE) SED_VERY_NEGATIVE;
      for (i=0;i<m->col;i++)
	{
	  x = (SED_MATTYPE) fabs ((double) m->mat[i][i]);
	  if ( (pivotcontrol[i] == TRUE)
	       && (x > max) )
	    { max = x; nr = i; }
	}
      pivotcontrol[nr] = FALSE;
      pivot = m->mat[nr][nr];
      if ( ((SED_MATTYPE) fabs ((double) pivot)) > nzero )
	{
	  ef++;
	  for (i=0;i<m->col;i++)
	    if (i!=nr)
	      {
		ai = m->mat[i][nr] / pivot;
		for (j=0;j<m->col;j++)
		  m->mat[i][j] -= ai*(m->mat[nr][j]);
		m->mat[i][nr] = ai;
	      }
	  for (j=0;j<m->col;j++)
	    m->mat[nr][j] /= (-pivot);
	  m->mat[nr][nr] = ((SED_MATTYPE) 1.0) / pivot;
	}
      else stop = TRUE;
    } 
  while ( (ef < m->col) && (stop == FALSE) );

  free (pivotcontrol);
  return (m->col - ef);
}

/*****************************************************************/

int SED_Tm_cholesky_decomposition (struct SED_longTmatrix *m)
/*
 *   labelled provisional in pascal version
 */
{
  SED_MATTYPE sum;
  int i,j,k;

  for (i=0;i<m->col;i++)
    for (j=i+1;j<m->col;j++)
      m->mat[i][j] = (SED_MATTYPE) 0.0;

  for (i=0;i<m->col;i++)
    {
      for (j=0;j<i;j++)
	{
	  sum = (SED_MATTYPE) 0.0;
	  for (k=0;k<j;k++) sum += (m->mat[i][k])*(m->mat[j][k]);
	  m->mat[i][j] = (m->mat[i][j] - sum) / m->mat[j][j];
	}
      sum = (SED_MATTYPE) 0.0;
      for (k=0;k<i;k++) sum += (m->mat[i][k])*(m->mat[i][k]);
      if (m->mat[i][i] - sum < 0) return -1;
      m->mat[i][i] = (SED_MATTYPE) sqrt (m->mat[i][i] - sum);
    }
  return 0;
}

/*****************************************************************/

int SED_Tm_choleskyinv (struct SED_longTmatrix *m)
{
  SED_MATTYPE sum, nzero;
  int i,j,k;
  int err=0;

#ifdef DEVEL_MAT_DEBUG
  if (m->row!=m->col)
    SED_Tm_error ("SED_Tm_choleskyinv", m, NULL, NULL, NULL, NULL);
#endif

  nzero = SED_Tm_trace (m);
  if (nzero!=0.0) nzero = nzero/((SED_MATTYPE)m->row) * NEARLYZERO;
  err=0;
  for (i=0;i<m->col;i++)
    {
      for (j=0;j<i;j++)
	{
	  sum = (SED_MATTYPE) 0.0;
	  for (k=0;k<j;k++) sum += (m->mat[i][k])*(m->mat[j][k]);
	  m->mat[i][j] = (m->mat[i][j] - sum) / m->mat[j][j];
	}
      sum = (SED_MATTYPE) 0.0;
      for (k=0;k<i;k++) sum += (m->mat[i][k])*(m->mat[i][k]);
      if (m->mat[i][i] - sum < 0) return -1;
      m->mat[i][i] = (SED_MATTYPE) sqrt(m->mat[i][i] - sum);
    }
  for (i=0;i<m->col;i++)
    if ( ((SED_MATTYPE) fabs ((double) m->mat[i][i])) < nzero)
      err++;
  if (err==0)
    {
      for (i=0;i<m->col;i++)
	{
	  for (j=0;j<i;j++)
	    {
	      sum = (SED_MATTYPE) 0.0;
	      for (k=j;k<i;k++) sum += (m->mat[i][k])*(m->mat[k][j]);
	      m->mat[i][j] = -(sum/m->mat[i][i]);
	    }
	  m->mat[i][i] = ((SED_MATTYPE) 1.0) / m->mat[i][i];
	}
      for (i=0;i<m->col;i++)
	for (j=0;j<=i;j++)
	  {
	    sum = (SED_MATTYPE) 0.0;
	    for (k=i;k<m->col;k++) sum += (m->mat[k][i])*(m->mat[k][j]);
	    m->mat[i][j] = m->mat[j][i] = sum;
	  }
    }
  return err;
}

/*****************************************************************/

SED_MATTYPE SED_Tm_sym_check (struct SED_longTmatrix *m)
{
  int i,j;
  SED_MATTYPE eps;

#ifdef DEVEL_MAT_DEBUG
  if (m->row!=m->col)
    SED_Tm_error ("SED_Tm_symcheck", m, NULL, NULL, NULL, NULL);
#endif
  eps = -1.0;
  for (i=0;i<m->col;i++)
    for (j=i;j<m->row;j++)
      if ( ((SED_MATTYPE) fabs ((double) (m->mat[i][j] - m->mat[j][i]))) > eps)
	eps = (SED_MATTYPE) fabs ((double) (m->mat[i][j] - m->mat[j][i]));
  return eps;
}

/*****************************************************************/

int SED_Tm_add_block (struct SED_longTmatrix *m,
		 struct SED_longTmatrix *block,
		 int row,
		 int col)
{
  int i,j;

#ifdef DEVEL_MAT_DEBUG
  if ( (row+block->row) > m->row
       || (col+block->col) > m->col )
    SED_Tm_error ("SED_Tm_add_block", m, block, NULL, &row, &col);
#endif

  for (i=0;i<block->row;i++)
    for (j=0;j<block->col;j++)
      m->mat[row+i][col+j] += block->mat[i][j];
  return 0;
}

/*****************************************************************/

int SED_Tm_pick_block (struct SED_longTmatrix *m,
		  struct SED_longTmatrix *block,
		  int row,
		  int col,
		  int drow,
		  int dcol)
{
  int i,j;

#ifdef DEVEL_MAT_DEBUG
  if ( (row==0) || (col==0) || (drow==0) || (dcol==0)
       || (m->row<row+drow) || (m->col<col+dcol)
       || (block->dim_row<drow) || (block->dim_col<dcol) )
    SED_Tm_error ("SED_Tm_pick_block", m, block, NULL, &row, &col);
#endif

  block->row = drow;
  block->col = dcol;
  for (i=0;i<drow;i++)
    for (j=0;j<dcol;j++)
      block->mat[i][j] = m->mat[row+i][col+j];
  return 0;
}

/*****************************************************************/

int SED_Tm_del_row (struct SED_longTmatrix *m,
	       int row)
{
  int i,j;

#ifdef DEVEL_MAT_DEBUG
  if ( (row>m->row) || (row<0) )
    SED_Tm_error ("SED_Tm_del_row", m, NULL, NULL, &row, NULL);
#endif

  /*
   *  shift up lower rows, if row is NOT lowest row
   *  thus upper limit is m->row - 1
   *  this bug fixed as in version 0.28
   */
  for (i=row;i<m->row-1;i++)
    for (j=0;j<m->col; j++)
      m->mat[i][j] = m->mat[i+1][j];
  m->row -= 1;
  return 0;
}

/*****************************************************************/


int SED_Tm_del_col (struct SED_longTmatrix *m,
	       int col)
{
  int i,j;

#ifdef DEVEL_MAT_DEBUG
  if ( (col>m->col) || (col<0) )
    SED_Tm_error ("SED_Tm_del_col", m, NULL, NULL, NULL, &col);
#endif

  /*
   *  shift left cols to the right, if col is NOT rightmost
   *  thus upper limit is m->col - 1
   *  this bug fixed as in version 0.28
   */
  for (i=col;i<m->col-1;i++)
    for (j=0;j<m->row; j++)
      m->mat[j][i] = m->mat[j][i+1];
  m->col -= 1;
  return 0;
}

/*******************************************************************/
/*
 *  weighted polynomial approximation using Least Mean Square methode
 *  f(x) = a + bx + cx**2 + dx**3 + ...
 */

int SED_Tm_polynom_approx (struct SED_longTmatrix *data,
		      int col_x,
		      int col_y,
		      int col_w,
		      const int col_app,
		      const int col_res,
		      const int row_beg,
		      const int row_end,
		      const int degree,
		      SED_MATTYPE x0,
		      SED_MATTYPE *dummy,
		      struct SED_longTmatrix *coeff,
		      struct SED_longTmatrix *qxx,
		      SED_MATTYPE *s0)
{
  int count;
  struct SED_longTmatrix *vec;
  SED_MATTYPE ya, v, vv, x, w;
  int i,j,k;

#ifdef DEVEL_MAT_DEBUG
  if ( (col_x>data->col-1) || (col_y>data->col-1) || (col_w>data->col-1) ||
       (col_app>data->col-1) || (col_res>data->col-1) ||
       (qxx->dim_row<degree) || (qxx->dim_col<degree) ||
       (coeff->dim_row<degree) || (row_end - row_beg + 1 < degree) ||
       (row_beg<0) || (row_beg>row_end) || (row_end>=data->row) )
    SED_Tm_error ("SED_Tm_polynom_approx", data, coeff, qxx, &col_x, &col_y);
#endif

  count = row_end - row_beg + 1;
  vec = SED_Tm_init (1, degree);
  qxx->col = qxx->row = coeff->row = degree;
  coeff->col = 1;
  SED_Tm_zero(qxx);
  SED_Tm_zero(coeff);
  vec->mat[0][0]=1;

  for (i=row_beg;i<=row_end;i++)
    if ( (dummy == NULL) || (data->mat[i][col_y] != *dummy) )
      {
	if (col_x < 0) x = ((SED_MATTYPE) i) - x0;
	else x = data->mat[i][col_x] - x0;
	if (col_w < 0) w = 1.0;
	else w = data->mat[i][col_w];
	for (j=1;j<degree;j++)
	  vec->mat[0][j] = vec->mat[0][j-1] * x;

	for (j=0;j<degree;j++)
	  {
	    for (k=0;k<degree;k++)
	      qxx->mat[j][k] += vec->mat[0][j]*vec->mat[0][k]*w;
	    coeff->mat[j][0] += data->mat[i][col_y]*vec->mat[0][j]*w;
	  }
      }

  SED_Tm_inv (qxx);
  SED_Tm_mult (qxx, FALSE, coeff, FALSE, coeff);

  vv = (SED_MATTYPE) 0.0;
  count = 0;
  for (i=row_beg;i<=row_end;i++)
    {
      if (col_x < 0) x = ((SED_MATTYPE) i) - x0;
      else x = data->mat[i][col_x] - x0;
      if (col_w < 0) w = 1.0;
      else w = data->mat[i][col_w];
      for (j=1;j<degree;j++)
	vec->mat[0][j] = vec->mat[0][j-1] * x;
      
      ya = (SED_MATTYPE) 0.0;
      for (j=0;j<degree;j++)
	ya += vec->mat[0][j]*coeff->mat[j][0];

      if ( (dummy == NULL) || (data->mat[i][col_y] != *dummy) )
	{
	  count++;
	  v = data->mat[i][col_y] - ya;
	  vv += v*v*w;
	  if (col_res >= 0) data->mat[i][col_res] = v;
	}
      else
	{ if (col_res >= 0) data->mat[i][col_res] = *dummy; }
      
      if (col_app >= 0) data->mat[i][col_app] = ya;
    }
  *s0 = (SED_MATTYPE) ((count - degree > 0) 
		  ? sqrt ((double) (vv / (count - degree))) 
		  : -999.0);
  SED_Tm_free (vec);
  return 0;
}

/*******************************************************************/
/*
 *  takes coefficients calculated by m_polynom_approx and calculates
 *  f(x).
 */
SED_MATTYPE SED_Tm_polynom_interpol (
	struct SED_longTmatrix *coeff,
        int num_coeff,
	const SED_MATTYPE dx)
{
	SED_MATTYPE res;
	int i;

#ifdef DEVEL_MAT_DEBUG
	if ((num_coeff-1)>coeff->row)
	  SED_Tm_error ("SED_Tm_polynom_interpol", coeff, NULL, NULL, &num_coeff, NULL);
#endif

	res = coeff->mat[num_coeff-1][0];

	for (i=num_coeff-2;i>=0;i--) {
		res = res * dx + coeff->mat[i][0];
	}
#ifdef LIB_DEBUG
	SED_Tm_write(stderr, coeff, 15, 8, 0, 0, 1);
	fprintf (stderr, "%15.8f\n", res );
#endif

	return res;
}

/*****************************************************************/
/*
 * separate matrix in symmetrical and asymmetrical part 
 * NULL pointer for matrix sym or asym indicates that this matrix is not calculated
 */
int SED_Tm_separate_sym(struct SED_longTmatrix *m, struct SED_longTmatrix *sym, struct SED_longTmatrix *asym)
{
  int i, j;

#ifdef DEVEL_MAT_DEBUG
  if (m->col!= m->row) 
    SED_Tm_error ("SED_Tm_separate_sym", m, NULL, NULL, &(m->row), &(m->col));
#endif

  for (i=0; i<m->col; i++)
    for (j=0; j<m->col; j++)
      {
	if (sym!=NULL)  sym->mat[i][j]  = (m->mat[i][j] + m->mat[j][i])/2;
	if (asym!=NULL) asym->mat[i][j] = (m->mat[i][j] - m->mat[j][i])/2;
      }
  return 0;
}


