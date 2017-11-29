/* pres.c */

/*  Brute force Kramer-Mesner matrix generation.
**
**  (c) 1997,1999 Kari J. Nurmela
**
** generate expression for a permutation while preserving the order
** of the bits
*/

#include <stdio.h>
#include <search.h>
#include <stdlib.h>
#include "pack.h"

#define pow_2(X) (((cword) 1) << (X))

struct moveRec {
  int indx, diff;
};

static int compareMoveRecs(const void* i, const void* j)
{
  return ((struct moveRec*)i)->diff - ((struct moveRec*)j)->diff;
}

int main(int argc, char *argv[])
{
  int from[CWORD_BITS], to[CWORD_BITS], i, n = argc-1;
  struct moveRec moves[CWORD_BITS+1];
  int start, end, ix;
  cword mask;

  /* read the arguments */

  if(n > CWORD_BITS) {
    fprintf(stderr, "Error: Too long permutation.\n");
    exit(1);
  }

  for(i = 1; i < argc; i++)
    from[i-1] = atoi(argv[i]) - 1;
  for(i = 0; i < n; i++)
    to[from[i]] = i;
  
  for(i = 0; i < n; i++) {
    moves[i].indx = i;
    moves[i].diff = to[i] - i;
  }

  qsort(moves, n, sizeof(struct moveRec), compareMoveRecs);
  moves[n].diff = CWORD_BITS;
  
  start = end = 0;
  while(start < n) {
    ix = start; mask = 0;
    do {
      mask |= pow_2(moves[ix++].indx);
    } while(moves[ix].diff == moves[ix-1].diff);
    end = ix - 1;
    if(moves[start].diff == 0)
      printf("((X) & %u)", mask);
    else if(moves[start].diff > 0)
      printf("((X) & %u) << %d", mask, moves[start].diff);
    else
      printf("((X) & %u) >> %d", mask, -moves[start].diff);
    start = ix;
    if(start < n)
      printf(" | \\\n");
  }
  return 0;
}
