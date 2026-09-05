#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define WITH_DEBUG
#define OS_LINUX

#include <common.h>

int main ()
{
	double x1=2.0;
	double x2=3.0;
	double p;
	double y1=10;
	double y2=20;
	double res;

	scanf ("%lf", &p);
	d1ipol_linear (x1,x2,y1,y2,p,&res);
	printf ("%12.8f\n", res);
}
