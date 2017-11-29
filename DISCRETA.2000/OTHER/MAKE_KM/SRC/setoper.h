/*  Brute force Kramer-Mesner matrix generation.
**
**  (c) 1997,1999 Kari J. Nurmela
**
**   setoper.h
**
**   This file contains declarations needed for subset ranking and unranking
**   plus making a complement of a given subset.
**
*/

#ifndef _setoper_h_
#define _setoper_h_

#include "bincoef.h"

/***** type definitions *****/

typedef unsigned short varietyType;
/* type of subset element numbers */
#define BIG_VARIETY ~((varietyType) 0)

/***** functions *****/

void getFirstSubset(varietyType *subset, int card);
int getNextSubset(varietyType *subset, int card, int v);

#endif
