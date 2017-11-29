/* utils.c */

/*  Brute force Kramer-Mesner matrix generation.
**
**  (c) 1997,1999 Kari J. Nurmela
*/

/* this file contains functions that DO NOT NEED perm.h */
/* this way the file pack.c is a bit smaller */

#include "pack.h"
#include "utils.h"
#include "setoper.h"
#include <stdio.h>
#include <stdlib.h>
#include <search.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <limits.h>

int len; /* codeword length */
int kPar;
int tPar;
int printTicks = 0;

int bits8[256];
cword power2[CWORD_BITS];

int no_search = 0;

char *KMfile = "";
int KMgroupOrder = 0;
int generatorCount = 0;

cword *code = 0;
long memoryLimit = LONG_MAX;

struct feasOrbitRec *feasOrbits = 0;
size_t feasOrbitCount = 0;
struct coverOrbitRec *coverOrbits = 0;
size_t coverOrbitCount = 0;
coverList_t *coverList = 0;
size_t coverListCount = 0;

int maxOrbitCount;

/* ======================== General utilities ======================= */

struct dynTableRec {
  size_t elemCount, elemSize, maxRef;
  void *memPtr;
};

static size_t dynUsedMemory = 0;

static void checkTableSize(struct dynTableRec *rec, size_t ref)
{
  if(ref >= rec->elemCount) {
    size_t newSize = 2 * rec->elemCount * rec->elemSize;
    if(newSize == 0)
      newSize = 4096; /* 4 K smallest allocation with checkTableSize */
    while(newSize < rec->elemSize * (ref + 1))
      newSize *= 2;
    if(printTicks) {printf("<%u>", newSize); fflush(stdout); }
    if((rec->memPtr = realloc(rec->memPtr, newSize)) == NULL) {
      fprintf(stderr, "Memory reallocation failed.\n");
      exit(1);
    }
    rec->elemCount = newSize / rec->elemSize;
  }
  if(ref > rec->maxRef) {
    dynUsedMemory += (ref - rec->maxRef) * rec->elemSize;
    rec->maxRef = ref;
  }
  if(dynUsedMemory > memoryLimit) {
    fprintf(stderr, "ML (%lu)\n", memoryLimit);
    exit(1);
  }
}

/* =========================== Tabular data ========================= */

void calculateBits8()
{
  int i, t;
  
  for(i = 0; i < 256; i++) {
    bits8[i] = 0; t = i;
    while(t) {
      if(t & 1) bits8[i]++;
      t >>= 1;
    }
  }
}

void calculatePower2()
{
  int i;

  power2[0] = 1;
  for(i = 1; i < CWORD_BITS; i++)
    power2[i] = 1 << i;
}

static cword minWordOnOrbit(cword cw)
{
  orbitSize_t orbitSize;
  cword minWord = cw;
  int i;

  orbitSize = generateOrbit(cw, orbit, MAX_ORBIT_LEN);
  if(generatorCount > 1)
    resetCtable(orbit, orbitSize);
  for(i = 0; i < orbitSize; i++)
    if(orbit[i] < minWord)
      minWord = orbit[i];
  return minWord;
}

void computeOrbitTables()
{
  struct baseOrbitRec *tmpOrbits = 0;
  int i, j;
  struct coverOrbitRec tmpCover, *coverPtr;
  int ones[CWORD_BITS+1], oneCount, card = tPar;
  varietyType subset[CWORD_BITS+1];
  cword coverWord;
  struct dynTableRec coverRec = {0, sizeof(coverList_t), 0, 0};

  free(coverList);
  coverList = 0;

  feasOrbitCount = makeOrbits(&tmpOrbits, kPar);
  feasOrbits = calloc(feasOrbitCount, sizeof(struct feasOrbitRec));
  if(!feasOrbits) {
    fprintf(stderr, "Memory allocation failed: feasOrbits.\n");
    exit(1);
  }
  for(i = 0; i < feasOrbitCount; i++) {
    feasOrbits[i].size = tmpOrbits[i].size;
    feasOrbits[i].firstWord = tmpOrbits[i].firstWord;
  }
  free(tmpOrbits); tmpOrbits = 0;
  if(printTicks) {
    printf("\nSorting %d feasible orbits...\n", feasOrbitCount);
    fflush(stdout);
  }
  qsort(feasOrbits, feasOrbitCount, sizeof(struct feasOrbitRec), 
	compareFeasOrbitRecs);
  if(printTicks) {
    printf("done\n");
    fflush(stdout);
    }

  coverOrbitCount = makeOrbits(&tmpOrbits, tPar);
  coverOrbits = calloc(coverOrbitCount, sizeof(struct coverOrbitRec));
  if(!coverOrbits) {
    fprintf(stderr, "Memory allocation failed: coverOrbits.\n");
    exit(1);
  }
  for(i = 0; i < coverOrbitCount; i++) {
    coverOrbits[i].size = tmpOrbits[i].size;
    coverOrbits[i].firstWord = tmpOrbits[i].firstWord;
  }
  free(tmpOrbits); tmpOrbits = 0;
  if(printTicks) {
    printf("\nSorting %d coverorbits...\n", coverOrbitCount);
    fflush(stdout);
  }
  qsort(coverOrbits, coverOrbitCount, sizeof(struct coverOrbitRec), 
	compareCoverOrbitRecs);
  if(printTicks) {
    printf("done\n");
    fflush(stdout);
  }

  /* overflows? */

  if((coverList_t) coverOrbitCount != coverOrbitCount) {
    fprintf(stderr, "Overflow: coverList_t\n");
    exit(1);
  }

  /* if(no_search)
     return; */

  if(generatorCount == 1)
    { free(ctable); ctable = NULL; }

  /* calculate covers */
  if(printTicks) { printf("Covers...\n"); fflush(stdout); }
  for(i = 0; i < coverOrbitCount; i++)
    coverOrbits[i].covered = 0;
  for(i = 0; i < feasOrbitCount; i++) {
    if(printTicks && i % 100 == 0) { printf(" %d ", i); fflush(stdout); }
    feasOrbits[i].coverCount = 0;
    for(j = 0, oneCount = 0; j < len; j++)
      if(power2[j] & feasOrbits[i].firstWord)
	ones[oneCount++] = j;
    if(oneCount != kPar) { fprintf(stderr, "error: kPar\n"); exit(1); }
    getFirstSubset(subset, card);
    do {
      for(j = 0, coverWord = 0; j < card; j++)
	coverWord |= power2[ones[subset[j]]];
      coverWord = minWordOnOrbit(coverWord);
      tmpCover.firstWord = coverWord;
      coverPtr = bsearch(&tmpCover, coverOrbits, coverOrbitCount,
			 sizeof(struct coverOrbitRec), compareCoverOrbitRecs);
      if(!coverPtr) { fprintf(stderr, "error: bsearch\n"); exit(1); }
      if(!coverPtr->covered) {
	/* coverPtr->covered = 1; allow multiple coverings */
	if(!feasOrbits[i].coverCount++)
	  feasOrbits[i].firstCover = coverListCount;
	checkTableSize(&coverRec, coverListCount + 1);
	coverList = (coverList_t *) coverRec.memPtr;
	coverList[coverListCount++] = coverPtr - coverOrbits;
      }
    } while(getNextSubset(subset, card, oneCount));
    /* for(j = 0; j < feasOrbits[i].coverCount; j++)
       coverOrbits[coverList[j + feasOrbits[i].firstCover]].covered = 0; */
    qsort(coverList + feasOrbits[i].firstCover, feasOrbits[i].coverCount,
	  sizeof(coverList_t), compareCoverList_ts);
  }
  if(printTicks) { printf("\ndone.\n"); fflush(stdout); }
}

void printKramerMesnerMatrix()
{
  covered_t *matrix;
  int i, j, co;
  FILE *fp;

  if(printTicks) printf("Writing Kramer-Mesner file...\n");
  matrix = calloc(feasOrbitCount*coverOrbitCount,
			     sizeof(covered_t));
  for(i = 0; i < feasOrbitCount; ++i)
    for(j = 0; j < feasOrbits[i].coverCount; ++j) {
      co = coverList[feasOrbits[i].firstCover+j];
      ++matrix[co+i*coverOrbitCount];
    }
  fp = fopen(KMfile, "w");
  if(!fp) {
    fprintf(stderr, "could not open file: %s\n", KMfile);
    exit(1);
  }
  fprintf(fp, "%% this file:    %s\n", KMfile);
  fprintf(fp, "%% order:        %d\n", KMgroupOrder);
  fprintf(fp, "%% degree:       %d\n", len);
  fprintf(fp, "%% # generators: %d\n", generatorCount);
  printGenerators(fp, len);
  fprintf(fp, "%% t, k:\n%% %d %d\n", tPar, kPar);
  fprintf(fp, "%% m, n:\n");
  fprintf(fp, "%d %d\n", coverOrbitCount, feasOrbitCount);
  for(j = 0; j < coverOrbitCount; ++j) {
    for(i = 0; i < feasOrbitCount; ++i)
      fprintf(fp, "%d ", matrix[j+i*coverOrbitCount] * feasOrbits[i].size /
	      coverOrbits[j].size);
    fprintf(fp, "\n");
  }
  fprintf(fp, "\nSTABILIZER-ORDER-K-SETS\n");
  for(i = 0; i < feasOrbitCount; ++i)
    fprintf(fp, "%d ", KMgroupOrder / feasOrbits[i].size);
  fprintf(fp, "\nSTABILIZER-ORDER-T-SETS\n");
  for(i = 0; i < coverOrbitCount; ++i)
    fprintf(fp, "%d ", KMgroupOrder / coverOrbits[i].size);
  fprintf(fp, "\n");
  fclose(fp);
  free(matrix);
  if(printTicks) printf("Matrix written\n");
}

/* ===================== ctable handling ====================== */

void setCtableBits(cword *orbit, orbitSize_t orbitSize)
{
  int i;
  cword *orbitPtr;

  for(i = 0, orbitPtr = orbit; i < orbitSize; i++, orbitPtr++)
    ctable[*orbitPtr / CWORD_BITS] |= 
      power2[*orbitPtr % CWORD_BITS];
}

void resetCtableBits(cword *orbit, orbitSize_t orbitSize)
{
  int k;
  cword *orbitPtr;

  for(k = 0, orbitPtr = orbit; k < orbitSize; k++, orbitPtr++)
    ctable[*orbitPtr / CWORD_BITS] &= ~power2[*orbitPtr % CWORD_BITS];
}

void resetCtable(cword *orbit, orbitSize_t orbitSize)
{
  int k;
  cword *orbitPtr;

  for(k = 0, orbitPtr = orbit; k < orbitSize; k++, orbitPtr++)
    ctable[*orbitPtr / CWORD_BITS] = 0;
}

void clearCtable()
{
  int i; cword *cptr;

  for(i = 0, cptr = ctable; i < ctableSize; i++)
    *cptr++ = 0;
}

/* ==================== Sorting functions ==================== */

int compareFeasOrbitRecs(const void* i, const void* j)
{
  return ((int) ((struct feasOrbitRec*)i)->firstWord) - 
    ((int) ((struct feasOrbitRec*)j)->firstWord);
}

int compareCoverOrbitRecs(const void* i, const void* j)
{
  return ((struct coverOrbitRec*)i)->firstWord - 
    ((struct coverOrbitRec*)j)->firstWord;
}

int compareCoverList_ts(const void* i, const void* j)
{
  return ((int) *((coverList_t *)i))  - ((int) *((coverList_t *)j));
}

/* ===================== orbit table generation ===================== */

#define TICK_MOD 10000

int makeOrbits(struct baseOrbitRec *orbits[], int orbWeight)
{

  /* generates orbits, 
     using weight orbWeight.
     Returns the number of generated orbits. */

  static int i, orbitSize;
  static cword firstw;
  static varietyType subset[CWORD_BITS+1], *ssptr;
  static int tickCounter = 0;
  int orbitCount = 0;
  struct dynTableRec orbitsRec = {0, sizeof(struct baseOrbitRec), 0, 0};

  free(*orbits);

  if(printTicks)
    printf("Number of words: %d\n", bincoef(len, orbWeight));

  fflush(stdout);
  getFirstSubset(subset, orbWeight);
  do {
    /* generate the first word of the orbit */
    firstw = 0;
    for(i = 0, ssptr = subset; i < orbWeight; i++, ssptr++)
      firstw |= power2[*ssptr];

    if(printTicks) {
      tickCounter++;
      while(tickCounter >= 2*TICK_MOD) {
	tickCounter -= 2*TICK_MOD;
	printf(".");
	fflush(stdout);
      }
    }

    /* generate the orbit, stop if already encountered */

    if(!(ctable[firstw / CWORD_BITS] & power2[firstw % CWORD_BITS])) {
      orbitSize = generateOrbitAndSetCtableBits(firstw, orbit, MAX_ORBIT_LEN);
      if(printTicks)
	tickCounter += orbitSize;
      checkTableSize(&orbitsRec, orbitCount + 1);
      *orbits = (struct baseOrbitRec*) orbitsRec.memPtr;
      (*orbits)[orbitCount].size = orbitSize;
      if((int) (*orbits)[orbitCount].size != (int) orbitSize) {
	fprintf(stderr, "Error: Too long orbit for orbits[].size.\n");
	exit(1);
      }
      (*orbits)[orbitCount].firstWord = firstw;
      
      orbitCount++;
    }

  } while(getNextSubset(subset, orbWeight, len));

  clearCtable();

  return orbitCount;
}

/* ===================== Argument parsing ===================== */

void parseArgs(int argc, char *argv[])
{
  int shift = 1;
  char *arg;

  len = atoi(argv[shift++]);
  kPar = atoi(argv[shift++]);
  tPar = atoi(argv[shift++]);
  while(shift < argc) {
    arg = argv[shift++];
    if(!strcmp("-pt", arg))
      printTicks = 1;
    else if(!strcmp("-ns", arg))
      no_search = 1;
    else if(!strcmp("-kmfile", arg)) {
      KMfile = argv[shift++];
      KMgroupOrder = atoi(argv[shift++]);
    }
    else {
      fprintf(stderr, "Error: Unknown argument: %s\n", arg);
      exit(1);
    }
  }
}

/* ============================ EOF ========================== */
