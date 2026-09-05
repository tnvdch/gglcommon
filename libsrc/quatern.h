/*
 *      FILE: quatern
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: quaternions
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
 *   AUTHORS: Cocard Marc <cocard@geod.ethz.ch>    
 *
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            18-09-1997  --  created
 */

/*
 *  now our own definitions
 */

#ifndef GGL_QUATER
#define GGL_QUATER

/* addition of two quaternions : result = q1 + q2 */
int q_add ( double q1[4], double q2[4], double result[4]);

/* multiplication of two quaternions: result = q1 * q2 */
int q_mult( double q1[4], double q2[4], double result[4]);

/* conjugate of quaternion 'q' is put in 'result' */
int q_conjugate( double q[4], double result[4]);

/* multipilcation of quaternion 'q' with scalar 'factor' put in 'result' */
int q_scalarmult ( double q[4], double factor, double result[4]);

/* returns the norm of the quaternion q */
double q_norm( double q[4]);

/* inverse of quaternion 'q' is put in 'result' */
int q_inverse( double q[4], double result[4]);

/* converts the rotation matrix 'rotmat' to the quaternion 'result' */
int rotmat_to_quat( double rotmat[3][3], double result[4]);

/* converts a rotation vector (must be a unit vector) and the associated
 * rotation angle (radian) to the corresponding quaternion 'result' */
int rotvec_to_quat( double rotvec[3], double angle, double result[4]);

/* converts the quaternion 'q' to the corresponding rotation matrix 'result' */
int quat_to_rotmat( double q[4], double result[3][3]);

/* rotates the vector 'x' using a quaternion 'q' to the rotated vector 'result' */
int q_rotate( double x[3], double q[4], double result[3]);

#endif

