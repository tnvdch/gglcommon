/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/chain.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: routines for keeping data in a chain structure
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: chain.c,v 1.10 2000/08/12 10:53:52 favey Exp favey $
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
 *   AUTHORS: Cocard Marc   <cocard@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            20-05-1998  --  created
 *
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/chain.c $
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.10  2000/08/12 10:53:52  favey
 * bugfix free_chain
 *
 * Revision 1.9  2000/07/07 14:51:49  cocard
 * *** empty log message ***
 *
 * Revision 1.8  1999/10/14 11:19:08  favey
 * bugfix in chain_init
 *
 * Revision 1.7  1999/09/08 09:32:31  favey
 * copyright info added
 *
 * Revision 1.6  1999/06/02 16:13:19  cocard
 * *** empty log message ***
 *
 * Revision 1.5  1998/08/25 16:18:28  cocard
 * c++ compatibility
 *
 * Revision 1.4  1998/08/05 17:24:24  cocard
 * small changes for c++ compilation
 *
 * Revision 1.3  1998/05/29 16:56:23  favey
 * *** empty log message ***
 *
 * Revision 1.2  1998/05/29 09:05:18  cocard
 * cosmetic
 *
 *
 */

/*
 * Standard Libaries
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "chain.h"

/*-----------------------------------------------------------------
  initialize a new chain
 */
struct Node * chain_init(int *count)
{
  *count=0;
  return chain_newelement(NULL);
}

/*-----------------------------------------------------------------
  allocates memory for a new element
*/
struct Node * chain_newelement(void * data)
{
  struct Node *chain;
   chain=(struct Node *) malloc(sizeof(struct Node));
#ifdef WITH_DEBUG
   if (chain==NULL)
     {
       fprintf(stderr, "GGL Common Library -- chain_init:\n");
       fprintf(stderr, "memory allocation error\n");
     }
#endif
   if (chain != NULL) { 
     memset(chain, 0, sizeof(struct Node));
     chain->data = data;
   }
   return chain;
}
/*-----------------------------------------------------------------
  free memory for an element
 */
void chain_freeelement(struct Node *p, int removedata)
{
  if (p==NULL) return;
  if (removedata && p->data!=NULL) free(p->data);
  free(p);
}

/*-----------------------------------------------------------------
  frees the  entire chain 
 */
void chain_free(struct Node * chain, int removedata, int  count)
{
  struct Node * p;
  
  chain=chain_head(chain);
  if (chain==NULL) return;
  
  if (chain->next!=NULL && chain->next->prev!=chain)
    chain->next->prev->next=NULL; /* this is for cyclic chains */

  while (chain!=NULL) { /* bugfix, we need to free also last element */
    p=chain;
    chain=chain->next;
    if (removedata && p->data!=NULL) free(p->data);
    free(p);
    count --;
#ifdef WITH_DEBUG
    if (count<0) {
      fprintf(stderr, "GGL Common Library -- chain_free:\n");
      fprintf(stderr, "count < 0\n");
    }
#endif
  }
}

/*-----------------------------------------------------------------
  inserts the member elememt in the chain (after chain)
 */
struct Node * chain_insert(struct Node * chain, struct Node * element, int *count)
{
  struct Node * p = chain->next;
  chain->next = element;
  element->prev =chain;
  element->next =p;
  if (p!=NULL) p->prev=element;
  (*count)++;
  return chain;
}

/*-----------------------------------------------------------------
  inserts the member element at the end of the chain
 */
struct Node * chain_append(struct Node * chain, struct Node * element, int *count)
{
   return chain_insert(chain_last(chain), element, count);
}

/*-----------------------------------------------------------------
  inserts the member element at the correct position using compare,
  return values of compare :
      positive  . if element1 > element2
      zero        if elememt1 = element2
      negative    if element1 < element2
 */
struct Node * chain_insertsorted( struct Node *chain, struct Node *element,
				  comparefunction compare, int *count)

{
   while (chain->next!=NULL && (*compare)(chain->next->data, element->data)<0) 
     chain=chain->next;
   while (chain->prev!=NULL && (*compare)(element->data, chain->data)<0)       
     chain=chain->prev;
   return chain_insert(chain, element, count);
}

/*-----------------------------------------------------------------
  remove actual element from chain
 */
struct Node * chain_remove(struct Node *chain, int *count)
{
  struct Node *p;
  if (chain==NULL || chain->prev==NULL)
    {
#ifdef WITH_DEBUG
       fprintf(stderr, "GGL Common Library -- chain_remove:\n");
       fprintf(stderr, "trying to remove head\n");
#endif
      return chain;
    }
  p=chain;
  chain=chain->prev;
  chain->next=p->next;
  chain->next->prev=chain;
  (*count)--;
  return chain;
}

/*-----------------------------------------------------------------
 returns the pointer to the head of the chain
 */
struct Node * chain_head(struct Node *chain)
{
  while (chain->prev!=NULL) chain=chain->prev;
  return chain;
}

/*-----------------------------------------------------------------
  returns the pointer to the last member in the chain
 */
struct Node * chain_last(struct Node *chain)
{
  while (chain->next!=NULL) chain=chain->next;
  return chain;
}

/*-----------------------------------------------------------------
  makes the chain cyclic by closing last with first element
  head is returned must be saved for freeing the chain
 */
struct Node * chain_cyclic(struct Node *chain)
{
  struct Node * head=chain_head(chain);
  struct Node * tail=chain_last(chain);
  if (head==tail) 
    {
#ifdef WITH_DEBUG
      fprintf(stderr, "GGL Common Library -- chain_cyclic:\n");
      fprintf(stderr, "trying to make zero chain cyclic\n");
#endif
      return head;
    }
  head->next->prev=tail;
  tail->next=head->next;
  return head;
}

/*-----------------------------------------------------------------
  makes the cyclic chain acylcic
 */
struct Node * chain_acyclic(struct Node *head)
{
#ifdef WITH_DEBUG
  if (head->prev!=NULL)
    {
      fprintf(stderr, "GGL Common Library -- chain_ucyclic:\n");
      fprintf(stderr, "head has previous element\n");
    }
  if (head->next==NULL)
    {
      fprintf(stderr, "GGL Common Library -- chain_ucyclic:\n");
      fprintf(stderr, "trying to make zero chain acyclic\n");
    }
#endif
  head->next->prev->next=NULL;  
  head->next->prev=head;
  return head;
}



/*-----------------------------------------------------------------
  returns the pointer to the next element of the chain
 */
struct Node * chain_next(struct Node *chain)
{
  if (chain==NULL) return NULL;
  return  chain->next;
}

/*-----------------------------------------------------------------
  returns the pointer to the previous element of the chain
 */
struct Node * chain_prev(struct Node *chain)
{
  if (chain==NULL) return NULL;
  return  chain->prev;
}


/*-----------------------------------------------------------------
  converts a chain to an array (with allocation of memory for array of pointers)
 */
void * chain_to_array(struct Node *chain, int count)
{
  int i;
  void  **array;

  array = (void **) malloc(sizeof(void *)*count);
#ifdef WITH_DEBUG
  if (array==NULL)
    {
      fprintf(stderr, "GGL Common Library -- chain_to_array:\n");
      fprintf(stderr, "memory allocation error\n");
    }
#endif
  chain=chain_head(chain);
  for (i=0;i<count;i++)
    {
      chain=chain->next;
      array[i]=chain->data;
    }  
  return array;
}

/*-----------------------------------------------------------------
  converts an array to a chain (with allocation of memory for chain)
 */
struct Node * array_to_chain(void **array, int count)
{
  struct Node *chain, *chain0;
  int idummy, i;

  chain0=chain=chain_init(&idummy);
  for (i=0; i<count; i++)
    chain=chain_append(chain, chain_newelement((array)[i]), &idummy);
  return chain0;
}

/*-----------------------------------------------------------------
  frees the array  
 */
void array_free(void **array, int removedata, int count)
{
  int i;
  if (removedata)
    for (i=0;i< count;i++) free(array[i]);    
  free(array);
}

/*-----------------------------------------------------------------
  sorts the elements from node1 to node2 using quicksort algorithm
 */
void chain_quicksort1(struct Node *node1,  struct Node *node2, comparefunction compare)
{
  struct Node *n1 = node1->prev;
  struct Node *n2 = node2;
  void *data = node2->data;
  void *dummy;
  
  for (;;)
    {
      do n1=n1->next;  while (n1!=n2 && compare(n1->data, data)<0);
      do n2=n2->prev;  while (n1!=n2 && compare(n2->data, data)>0);
      if (n1==n2 || n1->prev==n2) break;
      dummy=n1->data;  n1->data=n2->data;  n2->data=dummy;
    }
  dummy=n1->data;  n1->data=node2->data;  node2->data=dummy;
  if (node1!=n1 && node1!=n1->prev)  chain_quicksort1(node1, n1->prev, compare);
  if (node2!=n1 && n1->next!=node2)  chain_quicksort1(n1->next, node2, compare);
}  

/*-----------------------------------------------------------------
  sorts the entire chain and returns head of sorted chain 
*/
struct Node * chain_quicksort2(struct Node *chain, comparefunction compare)
{
  struct Node *head=chain_head(chain);
  chain_quicksort1(head->next, chain_last(chain), compare);
  return head;
}

/*-----------------------------------------------------------------
 */
