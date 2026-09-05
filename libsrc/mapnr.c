/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/mapnr.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: determine map number from swiss coordinates
 *
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: mapnr.c,v 1.3 2000/08/15 16:38:13 favey Exp $
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
 *   AUTHORS: Marc Troller (VTB Sommer 1999)
 *            Favey Etienne <favey@geod.ethz.ch>
 *
 *   HISTORY: 
 * $Log: mapnr.c,v $
 * Revision 1.3  2000/08/15 16:38:13  favey
 * mapras: bugfix and complete rewrite
 *
 * Revision 1.2  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.1  1999/06/30 15:34:56  favey
 * Initial revision
 *
 */
/* erstellt im SS99 waehrend des Vertiefungsblocks Flottenmanagement */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "errproc.h"
#include "mapnr.h"

int nrsheets25 (const double west,const double north,
		const double east,const double south,
		int *sheets, const int size)
/* 
 * Berechnet aus den Landeskartenblaetter der Eckpunkte einen 
 * array mit allen enthaltenden Kartenblaettern
 */
{
  int karte_1, karte_n, i, j, nr, anzy, anzx;
  double l1, l2, o1, o2;

  nr=0;
  karte_1=map_25(west, north);
  karte_n=map_25(east, south);
  if (karte_1 == karte_n) {
    sheets[0]=karte_1;
    return 1;
  } else {
    edge_25(karte_1, &l1, &o1);
    edge_25(karte_n, &l2, &o2);
    anzy=floor((l2-l1+1)/BREITE25)+1;
    anzx=floor((o1-o2+1)/HOEHE25)+1;
    for (j=1; j<=anzx; j++) for (i=1; i<=anzy; i++)
      if (nr>(size-1)) {
	puterr ("mapnr","nrsheets25");
	fprintf (stderr, "too many maps, buffer too small\n");
	return EXIT_FAILURE;
      } else {
	nr++;
	sheets[nr-1]=map_25(west+(i-1)*BREITE25,
			    north-(j-1)*HOEHE25);
      }  
    return nr;
  }
}

void edge_25 (int sheet_nr, double *links, double *oben)
/*Berechnet das linke obere Koordinatenpaar einer gegebenen Karte */
{
  int rechts, tief;

  rechts=sheet_nr-20*floor(sheet_nr/20.0);
  *links=BREITE25*rechts+480000;
  tief=(sheet_nr-1000)/20;
  *oben=-HOEHE25*tief+302000;
}


int map_25 (const double y, const double x) 
/*Berechnet aus den Koordinaten das entsprechende 1:25000er Kartenblatt*/ 
{
  double rechts, tief;

  rechts=floor((y-445000.0)/BREITE25)-2;
  tief=floor((338000.0-x)/HOEHE25)-3;
  return (1000+((int) rechts)+20*((int) tief)); 
}


int nrsheets100 (const double west,const double north,
		 const double east,const double south,
		 int *sheets, const int size)
/*
 * Berechnet aus den Landeskartenblaetter der Eckpunkte einen
 * array mit allen enthaltenden Kartenblaettern
 */
{
  int karte_1, karte_n, i, j, nr, anzy, anzx;
  double l1, l2, o1, o2;

  nr=0;
  karte_1=map_100(west, north);
  karte_n=map_100(east, south);
  if (karte_1 == karte_n) {
    sheets[0]=karte_1;
    return 1;
  } else {
    edge_100(karte_1, &l1, &o1);
    edge_100(karte_n, &l2, &o2);
    anzy=floor((l2-l1+1)/BREITE100)+1;
    anzx=floor((o1-o2+1)/HOEHE100)+1;
    for (j=1; j<=anzx; j++) for (i=1; i<=anzy; i++)
      if (nr>(size-1)) {
	puterr ("mapnr","nrsheets100");
	fprintf (stderr,
		 "too many maps, buffer too small\n");
	return EXIT_FAILURE;
      } else {
	nr++;
	sheets[nr-1]=map_100(west+(i-1)*BREITE100,
			     north-(j-1)*HOEHE100);
      }
    return nr;
  }
}

void edge_100 (int sheet_nr, double *links, double *oben)
/*Berechnet das linke obere Koordinatenpaar einer gegebenen Karte */
{
  int rechts, tief;

  rechts=sheet_nr-5*floor(sheet_nr/5.0);
  *links=BREITE100*rechts+480000;
  tief=(sheet_nr-25)/5;
  *oben=-HOEHE100*tief+302000;
}


int map_100 (const double y, const double x) 
/*Berechnet aus den Koordinaten das entsprechende 1:100000er Kartenblatt*/
{
  double rechts, tief;

  rechts=floor((y-480000.0)/BREITE100);
  tief=floor((302000.0-x)/HOEHE100);
  return (25+((int) rechts)+5*((int) tief));
}


