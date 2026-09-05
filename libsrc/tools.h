/*
 *      FILE: tools
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: 
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
 *   AUTHORS: Favey Etienne   <favey@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            31-08-1997  --  created
 *
 */

/*
 * Calculates an intersection of two arcs. The routine
 * gives the same result as intersection_circles from
 * unit intersect. Latter should be used. (note that
 * the return value of arc_intersect (= error) differ from
 * that of intersection_circles (=# of intersections).
 */

int arc_intersect (
  const double p1[3], /* x, y, radius of circle 1 */
  const double p2[3], /* x, y, radius of circle 2 */
  double na[2],       /* x, y of intersection 1 */
  double nb[2])       /* x, y of intersection 2 */
;

/*
 * initializes the first combination. All following
 * can be obtained by using combMofN_next
 */
int combMofN_first (int *comb, const int M)
;

/*
 * replaces comb by the next combination. returns 0 on success
 * and -1 after the last combination was reached. 
 * example :
 * combMofN_init(4,3)  comb = 0 1 2  (return value =  0)
 * combMofN_next(4,3)  comb = 0 1 3  (return value =  0)
 * combMofN_next(4,3)  comb = 0 2 3  (return value =  0)
 * combMofN_next(4,3)  comb = 1 2 3  (return value =  0)
 * combMofN_next(4,3)  comb = 1 2 3  (return value = -1 -> no new combination returned)
 */
int combMofN_next (int *comb, const int M, const int N)
;

/*
 * returns the number of combinations (M of N) as a double 
 */
double combMofN_count(const int M, const int N)
;

/*----------------------------------------------------------------------
 * routines for determining all permutation within an array
 * example for count = 3
 * 1:   0  1  2
 * 2:   1  0  2
 * 3:   0  2  1
 * 4:   2  0  1
 * 5.   1  2  0
 * 6:   2  1  0
 *  routines
 *      permutation_init  
 *      permutation_next 
 *      permutation_free
 *      permutation_ctrlprint
 *
 * ( for good understanding of how to use this routines see
 *   the code of permutation_ctrlprint )
 *---------------------------------------------------------------------- */


/*----------------------------------------------------------------------
 * initialize memory for the array <sequence>,  and for an additional
 * control array <ctrlarray> which is needed in permutation_next
 * n is the dimension of the array
 * fill the array with the first sequence
 */
 int permutation_init(int n, int **sequence, int **ctrlarray);

/*--------------------------------------------------------------------
 * permutation to the next sequence
 */
int permutation_next(int n, int * sequence, int * ctrlarray);

/*--------------------------------------------------------------------
 * returns the number of possible permutations
 * works up to n=12 then overflow
 */
long int permutation_count(int n);

/*--------------------------------------------------------------------
 * free the memory allocated by permutation_init
 */
int permutation_free(int *sequence, int *ctrlarray);

/*--------------------------------------------------------------------
 * routine for writing a complete set of all permutations 
 * of n integer in the range from 0..(n-1) to a file 
 * if addcounter == true the lines are numbered
 */
int permutation_ctrlprint(FILE *fp, int n, int addcounter);


