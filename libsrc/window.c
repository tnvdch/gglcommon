/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/window.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: file input and output functions including
 *            error handling
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: window.c,v 1.18 1999/09/08 09:12:18 favey Exp $
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
 *   AUTHORS: Favey Etienne   <favey@geod.ethz.ch>
 *
 *   HISTORY: 
 * $Log: window.c,v $
 * Revision 1.18  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.17  1999/06/30 13:21:50  favey
 * bugfixes with openzipped
 *
 * Revision 1.16  1999/06/30 12:17:30  favey
 * new function iszipped()
 *
 * Revision 1.15  1999/06/21 12:22:44  favey
 * bugfix in get_fname_wopath
 *
 * Revision 1.14  1999/06/21 09:47:52  favey
 * added file name to struct for easier logging
 *
 * Revision 1.13  1999/06/03 10:09:13  cocard
 * print_data_window: 2 format parameters added
 *
 * Revision 1.12  1999/04/08 13:32:52  favey
 * bugfix in time_index_data_window
 *
 * Revision 1.11  1998/12/15 17:00:11  favey
 * final behaviour of update_window (discussed with EF,MC)
 *
 * Revision 1.10  1998/12/14 17:33:12  favey
 * fixed several behaviour bugs in update_window
 *
 * Revision 1.9  1998/12/14 16:57:25  favey
 * no error on filling buffer in update_window!
 *
 * Revision 1.8  1998/10/30 14:22:07  favey
 * changed config.h to global.h
 *
 * Revision 1.7  1998/10/06 16:18:45  cocard
 * function time_index_data_window added
 *
 * Revision 1.6  1998/10/05 17:01:18  cocard
 * *** empty log message ***
 *
 * Revision 1.5  1998/10/02 18:47:17  cocard
 * major changes to readfunction
 *
 * Revision 1.4  1998/07/21 12:40:21  favey
 * added error function
 *
 * Revision 1.3  1998/05/07 14:39:30  favey
 * *** empty log message ***
 *
 * Revision 1.2  1998/05/04 07:02:14  favey
 * *** empty log message ***
 *
 * Revision 1.1  1998/04/29 10:28:58  favey
 * Initial revision
 *
 *
 */

/*
 *  standard C libraries
 */
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>

/*
 *  GGL Common Libs Headers
 */
#include "global.h"
#include "errproc.h"
#include "inout.h"
#include "matcalc.h"
#include "window.h"

/*------------------------------------------------------------------
 *  initialize windowing structure
 */
int init_data_window(
	const double half_window_size,
	const int buffer_size,
	const int numcols,
	const char *datafname,
	struct data_window *w,
	int (*read_new_element)(FILE *fpin, struct matrix *m, int row))

{
	int err;

	if ((w->m = m_init (buffer_size, numcols)) == NULL) return -3;
	if (buffer_size<2) return -5;
	if ((w->fpdata = openzipped (datafname, &(w->ispipe))) == NULL)
		return -2;
	w->fname = get_fname_wopath (datafname);

	if ( (err = (*read_new_element)(w->fpdata, w->m, 0)) != 1) 
		return (err==EOF)?EOF:-2;
	if ( (err = (*read_new_element)(w->fpdata, w->m, 1)) != 1) 
		return (err==EOF)?EOF:-2;
	w->m->row = 2;
	w->valid_count = 1;
	w->halfwinsize = half_window_size;
	w->bufsize = buffer_size;
	w->tail_time = w->head_time = w->mid_time = w->m->mat[0][0];

	return 0;
}

/*------------------------------------------------------------------
 *  update windowing structure
 */
int update_data_window(
	struct data_window *w,
	const double new_time,
	int (*read_new_element)(FILE *fpin, struct matrix *m, int row))
     
{
	double new_tail_time, new_head_time;
	int err, count, i, total=0;

	new_tail_time = new_time - w->halfwinsize;
	new_head_time = new_time + w->halfwinsize;

	i = w->valid_count;
	do {
		while (i>1 && w->m->mat[0][0] < new_tail_time)
		{
			m_pushup(w->m);
			i--;
		}

		count=0;
		while (w->m->mat[i][0] < new_head_time && i+1 < w->bufsize) 
		{
			i++;
			if ((err = (*read_new_element)(w->fpdata, w->m, i)) != 1) 
				return (err==EOF)?EOF:-2;
			count++;
		}
		total += count;
		w->m->row=i+1;
	} while (count > 0);

	w->valid_count=i;  
	w->tail_time = w->m->mat[0][0];
	w->head_time = w->m->mat[w->valid_count-1][0];
	w->mid_time  = (w->tail_time + w->head_time) / 2.0;

	if (new_time < w->tail_time) return -4;
	if (w->m->mat[i][0] < new_head_time) return -5;

	return total;
}

/*------------------------------------------------------------------
 *  eliminate windowing structure
 */
void eliminate_data_window (struct data_window *w)

{
	m_free (w->m);       w->m = NULL;
	if (w->ispipe) { pclose(w->fpdata); } else { fclose (w->fpdata); }
	w->fpdata = NULL;
	if (w->fname != NULL) free (w->fname);
	w->bufsize = w->valid_count = 0;
	w->mid_time = w->tail_time = w->head_time = w->halfwinsize = 0.0;
}

/*------------------------------------------------------------------
 * returns index i of row for  t(i) <= t < t(i+1)
 */
int time_index_data_window(struct data_window *w, double t)
{
	int i;

	for (i=0; i<w->valid_count-1; i++)
	{
		if (t>=w->m->mat[i][0] && t<w->m->mat[i+1][0]) return i;
	}
	return -1;
}

/*------------------------------------------------------------------
 *  debugging windowing structure
 */
void print_data_window (FILE *fpout, struct data_window *w, int fmt1, int fmt2)
{
  fprintf (fpout, "%14.7f %14.7f %14.7f %14.7f %5d %5d\n",
	   w->mid_time, w->tail_time, w->head_time,
	   w->halfwinsize, w->bufsize, w->valid_count);
  m_write (fpout, w->m, fmt1, fmt2, 0, 0, 1);
}

/*------------------------------------------------------------------
 *  print appropriate error message
 */
void window_err (const int errval, const char* fname)
{
	extern char *program_invocation_short_name;
	char *errmsg;
  
	fprintf (stderr, "%s: GGL Common Lib -- window\n",
		 program_invocation_short_name);
	switch (errval) 
	{
	case -1: errmsg = "end of file";    break;
	case -2: errmsg = "I/O error or incorrect matrix update by your function"; break;
	case -3: errmsg = "necessary memory could not be allocated";	break;
	case -4: errmsg = "time not in window"; break;
	case -5: errmsg = "buffer too small"; break;
	default: errmsg = "not a window error";
	}
	if (fname == NULL) {
		fprintf (stderr, "%s: %s\n", program_invocation_short_name,
		errmsg);
	} else {
		fprintf (stderr, "%s: %s: %s\n", program_invocation_short_name,
		fname, errmsg);
	}
}
