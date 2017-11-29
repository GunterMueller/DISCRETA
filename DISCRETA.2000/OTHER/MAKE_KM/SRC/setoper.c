/*  Brute force Kramer-Mesner matrix generation.
**
**  (c) 1997,1999 Kari J. Nurmela
**
**   setoper.c
**
**   This file contains functions needed for subset ranking and unranking
**   plus making a complement of a given subset.
**
*/


#include <stdio.h>
#include "pack.h"
#include "bincoef.h"
#include "setoper.h"


/*
** `getFirstSubset' gets the subset with rank 0 and the given cardinality.
**
*/

void getFirstSubset(varietyType *subset, int card)
{
  int i;

  for(i = 0; i < card; i++)
    *subset++ = i;
  *subset = BIG_VARIETY; /* sentinel */
}


/*
** `getNextSubset' gets the subset that has rank one bigger than the rank
** of the given subset. `v' is the number of varieties. Varieties range
** from 0 to v - 1. It returns 0, if this was the greatest possible rank
** for this type subset, otherwise the return value is 1.
*/

int getNextSubset(varietyType *subset, int card, int v)
{
  int i,j;

  if(subset[0] >= v - card)
    return 0;
  else {
    j = 0;
    while(subset[j + 1] <= subset[j] + 1)
      j++;
    subset[j]++;
    for(i = 0; i < j; i++)
      subset[i] = i;
    return 1;
  }
}
