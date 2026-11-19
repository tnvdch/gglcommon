/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/sed_matcalc.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: sed_matcalc.h,v 1.8 2000/05/12 09:22:53 favey Exp favey $
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
 *   AUTHORS: Cocard Marc     <cocard@geod.ethz.ch> (pascal version)
 *            Favey Etienne   <favey@geod.ethz.ch>  (c port)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *                  1987  --  created
 *            19-03-1997  --  pascal version to port
 *            07-04-1997  --  c port
 *            11-12-1997  --  procedure m_write debugged  
 *                            procedure m_writefmt added
 *                            procedure m_writefilled modified and completed (MC)
 *            01-07-1998  --  m_separate_sym added (MC)
 *            02-11-1998  --  introduced SED functionality
 *
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/sed_matcalc.h $
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.8  2000/05/12 09:22:53  favey
 * changed devel version to DEVEL_MAT_DEBUG
 *
 * Revision 1.7  2000/04/17 15:54:16  cocard
 * *** empty log message ***
 *
 * Revision 1.6  2000/02/25 12:57:58  ypeter
 * m_polynom_approx: modified function call, added possible weightedobs.
 * m_polynom_interpol: moved from winterpol to matcalc
 *
 * Revision 1.5  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.4  1999/07/29 13:22:28  cocard
 * $new function m_resize added (MC)
 *
 * Revision 1.3  1999/04/19 09:51:06  cocard
 * *** empty log message ***
 *
 * Revision 1.2  1999/04/07 11:45:32  cocard
 * read/write binary debugged
 *
 * Revision 1.1  1998/11/03 18:00:51  favey
 * Initial revision
 *
 *            
 *   CAVEATS: contains error messages exiting main()
 *            when compiled with DEVEL_MAT_DEBUG defined
 *            returning EXIT_FAILURE
 *
 */

#ifndef GGL_SED_TMATCALC
#define GGL_SED_TMATCALC

struct SED_longTmatrix {
	char *name;
	SED_MATTYPE **mat;	/* normal access via mat[][] */
	SED_MATTYPE *fast;	/* allows fast access to matrix */
	int row;
	int col;
	int dim_row;
	int dim_col;
};


/* writes some error on stderr */
void SED_Tm_error (
	const char *procname,
	struct SED_longTmatrix *m1, 
	struct SED_longTmatrix *m2,
	struct SED_longTmatrix *m3,
	int *row, 
	int *col);

/* initializes a matrix structure */
struct SED_longTmatrix * SED_Tm_init (int row, int col);

/* initializes the name of a matrix */
int SED_Tm_name (struct SED_longTmatrix *m, const char *newname);

/* initializes a matrix structure with name */
struct SED_longTmatrix * SED_Tm_init_name (int row, int col, char *newname);

/* frees the structure matrix */
void SED_Tm_free (struct SED_longTmatrix *m);


/* resizes the matrix without new memory allocation 
 * sets all values to zero if set_zero == TRUE
 */
int SED_Tm_resize (struct SED_longTmatrix *m, 
		   int row, int col, int set_zero);


/* fill one element at position (row, col) with value */
int SED_Tm_fill (struct SED_longTmatrix *m,
		 int row, int col, SED_MATTYPE value);

/* fill all elements of row  with value */
int SED_Tm_fill_row (struct SED_longTmatrix *m,
		     int row, SED_MATTYPE value);

/* fill all elements of column col  with value */
int SED_Tm_fill_col (struct SED_longTmatrix *m,
		     int col, SED_MATTYPE value);

/* copies upper part to lower part */
int SED_Tm_fillsym (struct SED_longTmatrix *m);

/* fills matrix completely with zero */
int SED_Tm_zero (struct SED_longTmatrix *m);

/* fills one row of matrix with zero */
int SED_Tm_zero_row (struct SED_longTmatrix *m,	int row);

/* makes the matrix to a unit matrix */
int SED_Tm_unit (struct SED_longTmatrix *m);

/* fills all elements of matrix with value */
int SED_Tm_set (struct SED_longTmatrix *m, SED_MATTYPE value);

/* fills all elements with random values in the range 0 - fac */
int SED_Tm_random (struct SED_longTmatrix *m, SED_MATTYPE fac);

/* deep copy from source to target */
int SED_Tm_copy (struct SED_longTmatrix *source,
		 struct SED_longTmatrix *target);

/* ????????? */
int SED_Tm_pushup (struct SED_longTmatrix *m);
int SED_Tm_pushdown (struct SED_longTmatrix *m);


int SED_Tm_max_col (
	struct SED_longTmatrix *m,
	const int col,
	SED_MATTYPE *max,
	int *row)
;
int SED_Tm_min_col (
	struct SED_longTmatrix *m,
	const int col,
	SED_MATTYPE *min,
	int *row)
;
int SED_Tm_max (
	struct SED_longTmatrix *m,
	SED_MATTYPE *max,
	int *row,
	int *col)
;
int SED_Tm_min (
	struct SED_longTmatrix *m,
	SED_MATTYPE *min,
	int *row,
	int *col)
;
int SED_Tm_absmax (
	struct SED_longTmatrix *m,
	SED_MATTYPE *max,
	int *row,
	int *col)
;
int SED_Tm_write (
	FILE *fpout,
	struct SED_longTmatrix *m,
	const int fmt1,
	const int fmt2,
	int n,
	const int tr,
	const int title)
;
int SED_Tm_writefmt (
	FILE *fpout,
	struct SED_longTmatrix *m,
	const char * fmt,
	int n,
	const int tr,
	const int title)
;
int SED_Tm_write_filled (
	FILE *fpout,
	struct SED_longTmatrix *m,
	const int tr,
	SED_MATTYPE nearlyzero)
;
int SED_Tm_writebin (
	FILE *fpout,
	struct SED_longTmatrix *m)
;
int SED_Tm_writebin2 (
	char * fname, 
	struct SED_longTmatrix *m)
;
int SED_Tm_readbin (
	FILE *fpin,
	struct SED_longTmatrix **m)
;
int SED_Tm_readbin2 (
	char * fname, 
	struct SED_longTmatrix **m)
;
int SED_Tm_readascii(
	char *fname, 
	struct SED_longTmatrix **m, 
	int row2read,
	int col2read)
;

int SED_Tm_add (
	struct SED_longTmatrix *m1,
	struct SED_longTmatrix *m2,
	struct SED_longTmatrix *m_res)
;
int SED_Tm_sub (
	struct SED_longTmatrix *m1,
	struct SED_longTmatrix *m2,
	struct SED_longTmatrix *m_res)
;
int SED_Tm_mult (
	struct SED_longTmatrix *m1,
	const int m1tr,
	struct SED_longTmatrix *m2,
	const int m2tr,
	struct SED_longTmatrix *m_res)
;
int SED_Tm_multx (
	struct SED_longTmatrix *m,
	SED_MATTYPE scalar,
	struct SED_longTmatrix *m_res)
;
SED_MATTYPE SED_Tm_trace (struct SED_longTmatrix *m)
;
int SED_Tm_trans (
	struct SED_longTmatrix *m,
	struct SED_longTmatrix *m_res)
;
int SED_Tm_atpa1 (
	struct SED_longTmatrix *a,
	struct SED_longTmatrix *p,
	struct SED_longTmatrix *atpa)
;
int SED_Tm_det (
        struct SED_longTmatrix *m, 
        SED_MATTYPE *det)
;

int SED_Tm_pivot (
	struct SED_longTmatrix *m,
	int row,
	int col)
;
int SED_Tm_inv (struct SED_longTmatrix *m)
;
int SED_Tm_cholesky_decomposition (struct SED_longTmatrix *m)
;
int SED_Tm_choleskyinv (struct SED_longTmatrix *m)
;
SED_MATTYPE SED_Tm_sym_check (struct SED_longTmatrix *m)
;
int SED_Tm_add_block (
	struct SED_longTmatrix *m,
	struct SED_longTmatrix *block,
	int row,
	int col)
;
int SED_Tm_pick_block (
	struct SED_longTmatrix *m,
	struct SED_longTmatrix *block,
	int row,
	int col,
	int drow,
	int dcol)
;
int SED_Tm_del_row (
	struct SED_longTmatrix *m,
	int row)
;
int SED_Tm_del_col (
	struct SED_longTmatrix *m,
	int col)
;
int SED_Tm_polynom_approx (
	struct SED_longTmatrix *data,
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
;

SED_MATTYPE SED_Tm_polynom_interpol (
	struct SED_longTmatrix *coeff,
        int num_coeff,
	const SED_MATTYPE dx)
;

int SED_Tm_separate_sym(
	struct SED_longTmatrix *m, 
	struct SED_longTmatrix *sym, 
	struct SED_longTmatrix *asym)
;
#endif
