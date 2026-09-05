/*
 *      FILE: chain
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: routines for keeping data in a chain structure
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
 *   AUTHORS: Cocard Marc     <cocard@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            20-05-1998  --  created
 */

/*
 * structure of a chain :
 *
 *  nil <- head <--> node1 <--> node2 ....  <--> last node  -> nil
 *          |          |          |                  |  
 *         nil       data1      data2              data
 *
 * for cyclic chains : node1 and last node are connected 
 * (but head->next still points to node1)
 */

/*
 * main definition of node
 */
#ifndef CHAIN
#define CHAIN

typedef int (*comparefunction)(void *, void *);

struct Node
{
  struct Node *prev;  /* pointer to previous node */
  struct Node *next;  /* pointer to next node */ 
  void *data;         /* void pointer to data */
};

struct Node * chain_init(int *count);           /* initialize a new chain */
struct Node * chain_newelement(void * data);    /* allocates memory for a new node */
/* free memory for an element */
void chain_freeelement (struct Node * p, int removedata);
void chain_free(struct Node *chain, int removedata, int count); /* frees the  entire chain */

/* inserts  the elememt in the chain (after chain) */
struct Node * chain_insert(struct Node *chain, struct Node *element, int *count);

/* inserts the element as last member in the chain */
struct Node * chain_append(struct Node *chain, struct Node *element, int *count);

/* inserts the  element at the correct position using compare,
 *  return values of compare :
 *      positive  . if element1 > element2
 *      zero        if elememt1 = element2
 *      negative    if element1 < element2
 */

struct Node * chain_insertsorted (struct Node *chain, struct Node *element, 
				 comparefunction compare, int *count);

/* remove actual element from chain */
struct Node * chain_remove (struct Node * chain, int *count);

struct Node * chain_head(struct Node *chain);  /* returns the pointer to the head of the chain */
struct Node * chain_last(struct Node *chain);  /* returns the pointer to the last node */
struct Node * chain_next(struct Node *chain);  /* next node */
struct Node * chain_prev(struct Node *chain);  /* previous node */

/*
  makes the chain cyclic by closing last with first element
  head is returned must be saved for freeing the chain
 */
struct Node * chain_cyclic(struct Node *chain);

struct Node * chain_acyclic(struct Node *head);   /*   makes the cyclic chain acylcic  */

/* converts a chain to an array (with allocation of memory for array of pointers) */
void * chain_to_array(struct Node * chain, int count);

/*  converts an array to a chain (with allocation of memory for chain)  */
struct Node * array_to_chain(void **array, int count);

void array_free(void **array, int removedata, int count);  /*  frees the array  */

/* sorts the elements from node1 to node2 using quicksort algorithm */
void chain_quicksort1(struct Node *node1,  struct Node *node2, comparefunction compare);

/* sorts the entire chain and returns head of sorted chain (uses chain_quicksort1)  */
struct Node * chain_quicksort2(struct Node *chain,  comparefunction compare);

#endif
