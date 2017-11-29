/* utils.h */

/*  Brute force Kramer-Mesner matrix generation.
**
**  (c) 1997,1999 Kari J. Nurmela
*/

#ifndef _utils_h_
#define _utils_h_

#include "pack.h"

/***** type definitions *****/

typedef short unsigned coverCount_t;
typedef short unsigned neighCount_t;
typedef short unsigned coveredByList_t;
typedef unsigned short covered_t;
typedef unsigned short coverList_t;

struct baseOrbitRec {
  orbitSize_t size;
  cword firstWord;
};

struct feasOrbitRec {
  orbitSize_t size;
  cword firstWord;
  coverCount_t coverCount;
  int firstCover;
};

struct coverOrbitRec {
  orbitSize_t size;
  cword firstWord;
  coverCount_t coveredByCount;
  int firstCoveredBy;
  covered_t covered;
};

extern struct feasOrbitRec *feasOrbits;
extern struct coverOrbitRec *coverOrbits;
extern coverList_t *coverList;

/***** functions *****/

extern void calculateBits8();
extern void calculatePower2();
extern void computeOrbitTables();
extern void printKramerMesnerMatrix();

extern void resetCtableBits(cword *orbit, orbitSize_t orbitSize);
extern void setCtableBits(cword *orbit, orbitSize_t orbitSize);
extern void resetCtable(cword *orbit, orbitSize_t orbitSize);
extern void clearCtable();

extern int makeOrbits(struct baseOrbitRec *orbits[], int);

extern int compareFeasOrbitRecs(const void*, const void*);
extern int compareCoverOrbitRecs(const void*, const void*);
extern int compareCoverList_ts(const void*, const void*);

extern void parseArgs(int, char *[]);

/***** global variables *****/

extern int len; /* codeword length */
extern int kPar;
extern int tPar;
extern int printTicks;
extern int solCount;
extern int solModCount;
extern int bits8[];
extern cword power2[];
extern struct searchLimitRec searchLimits[];
extern int searchLimitCount;
extern int no_search, automatic, automax, autoMaxCombs;
extern int maxtime;
extern float tabuLenFact[];
extern int tabuLenFactCount;
extern int lastTabuLen;
extern int maxExhaustive;
extern FILE *resultFp;
extern int PRNGseed, PRNGseedGiven;
extern char *KMfile;
extern int printDistanceDistribution;
extern int printDistanceDistribution;
extern int multiStartCount;
extern int stopIfFound;
extern int continuedSearch;
extern char *infoString;
extern int tabuStrategy;
extern int windowLength;
extern int tabuGoal;
extern int limit01;
extern int rrwGoal;
extern int extendCode;
extern int extendLocalCount;
extern int extendKludge;
extern int quitEarly;
extern int printCycLimit;
extern int printBounds;
extern int printLatex;
extern int generatorCount;
extern int hits[];
extern cword *code;
extern int words;
extern int maxOrbitCount;
extern int stepCounter;

#define MAX_PROB_COUNT 100000


#endif
