/* pack.c */

/*  Brute force Kramer-Mesner matrix generation.
**
**  (c) 1997,1999 Kari J. Nurmela
*/

#include "pack.h"
#include "bincoef.h"
#include "setoper.h"
#include "perm.h"
#include "utils.h"
#include <stdlib.h>

#define CHAR_BITS 8

/* permutation generation globals */

static cword permtable[PERM_GEN_COUNT];
cword tmpcword;

/* orbit generation globals */

cword *ctable;
int ctableSize;
#if PERM_GEN_COUNT > 1
#define GENSTACK_SIZE MAX_ORBIT_LEN
static cword genStack[GENSTACK_SIZE];
int genStackSize;
#endif

/* misc. definitions */

cword orbit[MAX_ORBIT_LEN];

/*****************************************************************/

int generateOrbit(cword firstw, cword* orbit, int maxOrbitLen)
{

  /* generates an orbit. If multiple generators, leaves ctable bits on.
     If single generator, does not touch ctable */

#if PERM_GEN_COUNT == 1
  /* do not use ctable */
  cword cw, *orbitPtr, *lastPos;

  cw = orbit[0] = firstw;
  orbitPtr = orbit;
  lastPos = orbit + maxOrbitLen - 1;
  while(1) {
    PERMUTE(cw);
    cw = *permtable;
    if(cw == firstw)
      break;
    *++orbitPtr = cw;
    if(orbitPtr >= lastPos) {
      fprintf(stderr, "Too long orbit (see MAX_ORBIT_LEN).\n");
      exit(1);
    }
  }
  return (orbitPtr-orbit)+1;
#else
  cword cw, *orbitPtr, *permPtr, *lastPos;
  int i;
  /* assume that ctable entries for this orbit are zero,
     sets the entries to 1 while performing a depth-first search */
  genStack[0] = orbit[0] = firstw;
  genStackSize = 1;
  orbitPtr = orbit + 1;
  lastPos = orbit + maxOrbitLen - 1;
  ctable[firstw / CWORD_BITS] |= power2[firstw % CWORD_BITS];

  while(genStackSize) { /* depth-first search */
    cw = genStack[--genStackSize];
    PERMUTE(cw);
    for(i = 0, permPtr = permtable; i < PERM_GEN_COUNT; i++, permPtr++) {
      cw = *permPtr;
      if(!(ctable[cw / CWORD_BITS] & power2[cw % CWORD_BITS])) {
	genStack[genStackSize++] = *orbitPtr++ = cw;
	ctable[cw / CWORD_BITS] |= power2[cw % CWORD_BITS];
	if(orbitPtr > lastPos) {
	  fprintf(stderr, "Too long orbit (see MAX_ORBIT_LEN).\n");
	  exit(1);
	}
      }
    }
  }

  return orbitPtr-orbit;
#endif
}

int generateOrbitAndSetCtableBits(cword cw, cword* orbit, 
				  int maxOrbitLen)
{

  /* generates an orbit and sets the corresponding bits in ctable */

#if PERM_GEN_COUNT == 1
  int ol;
  ol = generateOrbit(cw, orbit, maxOrbitLen);
  setCtableBits(orbit, ol);
  return ol;
#else
  return generateOrbit(cw, orbit, maxOrbitLen);
#endif
}

void printGenerators(FILE *fp, int printCycLimit)
{
  int to[CWORD_BITS], i, j, pt, tmp;
  cword cw;

  for(j = 0; j < PERM_GEN_COUNT; j++) {
    for(i = 0; i < printCycLimit; i++)
      to[i] = 0;
    for(i = 0; i < printCycLimit; i++) {
      cw = power2[i];
      PERMUTE(cw);
      while(permtable[j] != 1) {
	permtable[j] >>= 1;
	to[i]++;
      }
    }
    fprintf(fp, "%% ");
    for(i = 0; i < printCycLimit; i++)
      if(to[i] != -1 && to[i] != i) {
	pt = to[i];
	fprintf(fp, "(%d", i+1);
	tmp = to[i]; to[i] = -1;
	while(pt != i) {
	  if(pt >= printCycLimit) {
	    fprintf(stderr, "printCycLimit incorrect\n");
	    exit(1);
	  }
	  fprintf(fp, ",%d", pt+1);
	  tmp = to[pt]; to[pt] = -1;
	  pt = tmp;
	}
	fprintf(fp, ")");
      }
    fprintf(fp, "\n");
  }
}

int main(int argc, char *argv[])
{
  calculateBinCoefs();
  calculateBits8();
  calculatePower2();

  parseArgs(argc, argv);

  generatorCount = PERM_GEN_COUNT;

  ctableSize = (1 << len) / CWORD_BITS;
  ctable = (cword*) calloc(ctableSize, sizeof(cword));
  if(!ctable) {
    fprintf(stderr, "Error: Could not perform calloc\n");
    exit(1);
  }

  computeOrbitTables();

  if(strlen(KMfile) > 0)
    printKramerMesnerMatrix();

  free(ctable); ctable = 0;

  return 0;
}
      
