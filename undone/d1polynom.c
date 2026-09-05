/*
 *      FILE: $Source: /mnt/molnya/linux/gglutils/src/RCS/common/interpol/d1polynom.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: Polynominal one dimensional interpolation
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: d1polynom.c,v 1.2 1998/04/28 16:02:06 favey Exp $
 *
 *   AUTHORS: Favey Etienne   (c version)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            28-04-1998  --  created
 * $Log: d1polynom.c,v $
 * Revision 1.2  1998/04/28 16:02:06  favey
 * *** empty log message ***
 *
 *
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
 *  our special headers
 */
#include "config.h"
#include "errproc.h"
#include "d1polynom.h"



