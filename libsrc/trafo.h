/**
 **        FILE: trafo.h
 **
 ** DESCRIPTION: Functions to transform coordinates
 **              from the bessel ellipsoid to the
 **              Swiss military coordinates and the
 **              vice versa.
 **              Constants and ellipsoid data for calculations
 **              with the swiss projection system.
 **
 **              This file contains the declarations for trafo.c
 **
 **   COPYRIGHT: Copyright (C) 1996 by Etienne Favey
 **
 **  DISCLAIMER:   This program is free software; you can redistribute
 **              it and/or modify it under the terms of the GNU General
 **              Public License as published by the Free Software
 **              Foundation; either version 2 of the License, or
 **              (at your option) any later version.
 **                This program is distributed in the hope that it will
 **              be useful, but WITHOUT ANY WARRANTY; without even the
 **              implied warranty of MERCHANTABILITY or FITNESS FOR A
 **              PARTICULAR PURPOSE. See the GNU General Public License
 **              for more details.
 **                You should have received a copy of the GNU General
 **              Public License along with this program; if not, write
 **              to the Free Software Foundation, Inc., 675 Mass Ave,
 **              Cambridge, MA 02139, USA.
 **
 **      AUTHOR: Etienne Favey <eti@favey.lugs.ch>
 **              Zurich, Switzerland
 **
 **     VERSION: 4
 **
 **     CHANGES: dd-mm-yyyy -- description
 **              03-12-1996 -- first release (Version 1.05)
 **              18-12-1996 -- bugfix xyz2blh iteration criterion
 **              10-04-1997 -- adjusted to use with GGL coordtr.c
 **                            taken from public version 1.08,
 **                            renamed blh to llh, introduced
 **                            struct ellipsoid and vectors
 **              07-05-1997 -- added admsstr(), ellsupport, det_ell()
 **              14-01-1998 -- added struct frameconv, select_ellipsoid(),
 **                            select_frameconv()
 **              20-01-1998 -- split swiss.c into swiss.c and trafo.c
 **
 **/
#ifndef TRAFO
#define TRAFO


/*
 *  Frame Conversion Structure to simplify datum shifts
 */
struct frameread {
	int fwd;	/* sign for forward or reverse conversion */
	double ref_epoch;
	double dx[3];	/* translation in [m] */
	double scale;	/* scale - 1 in [ppm] */
	double rot[3];	/* small rotation elements [rad] */
	double ddx[3];
	double dscale;	/* rates of the above per year */
	double drot[3];
} ;


struct frameconv {
	int fwd;
	double scale;
	double rm[3][3];
	double dx[3];
} ;


/*
 * Ellipsoid patameters for several frequently used frames
 */
struct ellipsoid {
  double a;
  double e2;
} ;

/*
 * conversion of cartesian coordinates to spherical
 */
void cartesian_to_spherical(
			    double xyz[3], 
			    double *lat, 
			    double *lon, 
			    double *radius);

/*
 * conversion of spherical to cartesian coordinates
 */
void spherical_to_cartesian(
			    double lat, 
			    double lon, 
			    double radius, 
			    double xyz[3]);

/*
 *  geocentric to ellipsoidal
 *  parameters have to be passed and are returned in rad
 */ 
int xyz_to_llh(
	       const struct ellipsoid ell,
	       const double xyz[3],
	       double *lat,
	       double *lon,
	       double *h);


/*
 *  ellipsoidal to geocentric
 *  parameters have to be passed and are returned in rad
 */
int llh_to_xyz(
	       const struct ellipsoid ell,
	       const double lat,
	       const double lon,
	       const double h,
	       double xyz[3]);


/*
 *  select ellipsoid from configuration file
 *  returns pointer to struct or NULL upon failure
 *  free struct with free(3) after finishing using it
 */
struct ellipsoid * select_ellipsoid (
	const char *ellstr)
;


/*
 *  transform arbitrary frame given source epoch
 */
int frame_conv (
	struct frameconv *fcp,
	double source[3],
	double target[3])
;


/*
 *  select frame conversion parameters from file
 *  returns pointer to struct or NULL upon failure
 *  free struct with free(3) after finishing using it
 */
struct frameconv * select_frameconv (
	const char *source,
	const double inepoch,
	const char *target,
	const double outepoch)
;

/*
 *  select form types for input
 *  see trafo config file for more information
 */
int select_formtype (
	const char *formtype,
	int ft_type[3],
	int ft_cols[3])
;

/*
 *  select geoid from trafo config file
 *  returns pointer to struct or NULL upon failure
 *  free struct with free(3) after finishing using it
 */
struct Grid * select_geoid (
	const char *geoidstr,
	int g_ftyp[3])
;

#endif

