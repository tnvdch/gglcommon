/*
 *      FILE: $Source: /mnt/moho/linux/gglutils/src/RCS/common/netcdf_grid.c $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: grid operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: netcdf_grid.c 0.3 1997/12/24 12:11:44 favey Exp $
 *
 *   AUTHORS: Favey Etienne   (c version)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            23-11-1997  --  created
 * $Log: netcdf_grid.c $
 * Revision 0.3  1997/12/24 12:11:44  favey
 * followup of V. 0.29
 *
 *
 * DOCUMENTATION: pipe this file into
 *
 *            grep " \*% " | cut -c 5- | nroff -man | less -is
 *
 *            to view the manpages describing the functions
 *            (That's how I do it under Linux)
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
 *  the netCDF library: http://www.unidata.ucar.edu/packages/netcdf/index.html
 */
#include <netcdf/netcdf.h>

/*
 *  our special headers
 */
#include <errproc.h>
#include <mat3.h>
#include <matcalc.h>
#include <inout.h>
#include <gridop.h>




int g_write_netcdf (
	const char *ncfname,
	struct GGL_grid *g)
{
	int status;
	int ncid;
	int grxid, gryid; /* dimensions */
	int grdimids[2];
	int grdatid;
	int nxid, nyid, xminid, xmaxid, yminid, ymaxid, zminid, zmaxid;
	static size_t start[] = {0, 0}; /* start at first value */
	static size_t count[2];
	count[0] = g->nx;
	count[1] = g->ny;
	

	status = nc_create (ncfname, NC_NOCLOBBER, &ncid);  /* open dataset */
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_def_dim(ncid, "dgrid_nx", g->nx, &grxid);
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}
	status = nc_def_dim(ncid, "dgrid_ny", g->ny, &gryid);
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}
	grdimids[0] = grxid;
	grdimids[1] = gryid;

	status = nc_def_var (ncid, "grid_nx", NC_INT, 0, 0, &nxid);
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_def_var (ncid, "grid_ny", NC_INT, 0, 0, &nyid);
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_def_var (ncid, "grid_xmin", NC_DOUBLE, 0, 0, &xminid);
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_def_var (ncid, "grid_xmax", NC_DOUBLE, 0, 0, &xmaxid);
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_def_var (ncid, "grid_ymin", NC_DOUBLE, 0, 0, &yminid);
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_def_var (ncid, "grid_ymax", NC_DOUBLE, 0, 0, &ymaxid);
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_def_var (ncid, "grid_zmin", NC_DOUBLE, 0, 0, &zminid);
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_def_var (ncid, "grid_zmax", NC_DOUBLE, 0, 0, &zmaxid);
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_def_var (ncid, "grid", NC_DOUBLE, 2, grdimids, &grdatid);
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_enddef(ncid);      /* put in define mode */
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_put_var_int (ncid, nxid, &(g->nx));
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_put_var_int (ncid, nyid, &(g->nx));
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_put_var_double (ncid, xminid, &(g->xmin));
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_put_var_double (ncid, xmaxid, &(g->xmax));
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_put_var_double (ncid, yminid, &(g->ymin));
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_put_var_double (ncid, ymaxid, &(g->ymax));
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_put_var_double (ncid, zminid, &(g->zmin));
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_put_var_double (ncid, zmaxid, &(g->zmax));
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	status = nc_put_vara_double (ncid, grdatid, start, count,
				     g->data->fast);
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}


	status = nc_close(ncid);       /* close netCDF dataset */
	if (status != NC_NOERR) {
#ifdef WITH_ERRORS
		fprintf(stderr, "%s\n", nc_strerror(status));
#endif
		return GGLERR_INOUT;
	}

	return 0;
}
