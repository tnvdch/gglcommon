#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <common.h>

#define RHO (180.0/M_PI)

static void rotdx (const char *label, const double v[3], 
		   const double yaw, const double roll)
{
	double tmp[3];
	rotate (v, 3, yaw, tmp);
	rotate (tmp, 2, roll, tmp);
/*  	rotate (tmp, 2, M_PI, tmp); */
/*  	rotate (tmp, 3, M_PI/2.0, tmp); */
	printf ("%s: %10.6f %10.6f %10.6f\n", label, tmp[0], tmp[1], tmp[2]);
}

static void calca (double left[3], double rght[3], double *yaw, double *roll)
{
	double dlr[3], tmp[3];
	v3_sub(rght, left, dlr);
	v3_copy(dlr, tmp); tmp[2] = 0.0;
	*roll = asin(-dlr[2]/v3_length(dlr));
	*yaw = acos (dlr[0]/v3_length(tmp));
	printf ("roll %10.8f yaw %10.8f\n", *roll*RHO, 360.0-*yaw*RHO);
	rotdx ("dlr ", dlr, *yaw, *roll);
}

void main ()
{
	double rght[3] = { 0.397,  4.448, 0.317};
	double left[3] = {-4.182,  1.352, 0.408};
	double frnt[3] = {-1.949,  3.416, 0.123};
	double tail[3] = { 2.843, -4.538, 2.824};
	double dlr[3];
	double roll, yaw;

	calca(left, rght, &yaw, &roll);
	rotdx("frnt", frnt, yaw, roll);
	rotdx("tail", tail, yaw, roll);
	rotdx("left", left, yaw, roll);
	rotdx("rght", rght, yaw, roll);
	v3_sub(rght, left, dlr);
	v3_norm(tail, tail); v3_norm(dlr,dlr);
	printf ("%10.5f\n", acos(v3_scalar(tail, dlr))*RHO);
}
