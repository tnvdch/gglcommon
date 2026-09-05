/*
 *  reads a surf grid format and writes data to 
 *  standard output which can be piped to xyz2grd
 *  thus creating a gmt grid file
 *  see also scriptfile    surf2gmt
 *
 *  MC 16.04.98
 */

#define OS_LINUX

#include <stdlib.h>
#include <stdio.h>
#include "common.h"


int main (int argc, char **argv)
{
  struct Grid * data;
  int i, j;

  if ( (data = surf_readgrid(argv[1])) ==NULL)  
    {
      fprintf(stderr,"Error in g_read  program aborted\n");
      return -1;
    }
  fprintf(stdout,"#COMMENT %s\n",argv[1]);
  fprintf(stdout,"#MINLAT  %f\n",data->ymin);
  fprintf(stdout,"#MAXLAT  %f\n",data->ymax);
  fprintf(stdout,"#MINLON  %f\n",data->xmin);
  fprintf(stdout,"#MAXLON  %f\n",data->xmax);
  fprintf(stdout,"#DLAT    %f\n",data->dy);
  fprintf(stdout,"#DLON    %f\n",data->dx);

  for (i=data->ny-1; i>=0; i--)
    {
       for (j=0; j<data->nx; j++)
	   fprintf(stdout,"%f\n", g_getvalue (data, j, i));
    }
  return EXIT_SUCCESS;
}
