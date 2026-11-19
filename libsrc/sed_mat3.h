/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/sed_mat3.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 * $Id: sed_mat3.h,v 1.4 2001/02/21 13:28:19 cocard Exp favey $
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *                  1987  --  created
 *            19-03-1997  --  pascal version to port
 *            07-04-1997  --  c port
 *            14-05-1997  --  fixed boundary bug, thanks Mario!
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/sed_mat3.h $
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.4  2001/02/21 13:28:19  cocard
 * *** empty log message ***
 *
 * Revision 1.3  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.2  1999/02/01 17:04:08  favey
 * initial update
 *
 * Revision 1.1  1999/02/01 17:02:54  favey
 * Initial revision
 *           
 *
 */

/*
 *  now our own definitions
 */
#ifndef GGL_SED_TMAT3
#define GGL_SED_TMAT3

int SED_Tm3_add            ( SED_MATTYPE a[3][3], SED_MATTYPE b[3][3], SED_MATTYPE c[3][3]);
int SED_Tm3_sub            ( SED_MATTYPE a[3][3], SED_MATTYPE b[3][3], SED_MATTYPE c[3][3]);
int SED_Tm3_unit           ( SED_MATTYPE a[3][3]);
int SED_Tm3_copy           ( SED_MATTYPE a[3][3], SED_MATTYPE c[3][3]);
int SED_Tm3_trans          ( SED_MATTYPE a[3][3], SED_MATTYPE c[3][3]);
int SED_Tm3_inv            ( SED_MATTYPE a[3][3]);
void SED_Tm3_write         ( FILE *fpout, SED_MATTYPE a[3][3]);

int SED_Tv3_add            ( SED_MATTYPE a[3], SED_MATTYPE b[3], SED_MATTYPE c[3]);
int SED_Tv3_sub            ( SED_MATTYPE a[3], SED_MATTYPE b[3], SED_MATTYPE c[3]);
int SED_Tv3_copy           ( SED_MATTYPE a[3], SED_MATTYPE c[3]);
int SED_Tv3_multx          ( SED_MATTYPE a[3], SED_MATTYPE x, SED_MATTYPE c[3]);
int SED_Tv3_norm           ( SED_MATTYPE a[3], SED_MATTYPE c[3]);
SED_MATTYPE SED_Tv3_length ( SED_MATTYPE a[3]);
SED_MATTYPE SED_Tv3_scalar ( SED_MATTYPE a[3], SED_MATTYPE b[3]);
int SED_Tv3_cross          ( SED_MATTYPE a[3], SED_MATTYPE b[3], SED_MATTYPE c[3]);

int SED_Tm3v3_mult         ( SED_MATTYPE a[3][3], SED_MATTYPE b[3], SED_MATTYPE c[3]);
int SED_Tm3m3_mult         ( SED_MATTYPE a[3][3], SED_MATTYPE b[3][3], SED_MATTYPE c[3][3]);

int SED_Trot_getmat        ( int axis, SED_MATTYPE angle, SED_MATTYPE rotmat[3][3]);
int SED_Trotate            ( SED_MATTYPE x1[3], int axis, SED_MATTYPE angle, SED_MATTYPE x2[3]);
int SED_Trotvec_to_mat     ( SED_MATTYPE rotvec[3], SED_MATTYPE angle,SED_MATTYPE rotmat[3][3]);
int SED_Trotmat_to_vec     ( SED_MATTYPE rotmat[3][3], SED_MATTYPE rotvec[3], SED_MATTYPE *angle);

#endif

