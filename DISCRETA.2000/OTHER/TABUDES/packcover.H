// packcover.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

// this file contains simple declarations and definitions that do not
// depend on declarations in other files

#ifndef PACKCOVER_packcover_h
#define PACKCOVER_packcover_h

#include <stdlib.h>
#include <time.h>

#ifdef PACKCOVER_SMALL_PATCHES
typedef short patchType;
#else
typedef int patchType;
#endif

#ifdef PACKCOVER_SMALL_POINTS
typedef short pointType;
#else
typedef int pointType;
#endif

//typedef bool distType;
typedef unsigned char distType;

typedef short lambdaType;

typedef int volumeType;
// volume always non-negative, but volume change can assume negative
// values, too

typedef int penaltyType;
// penalty always non-negative

typedef double impactType;
// impact of a move

#define DEFAULT_MAX_MEM 16777216L

/***** Random numbers *****/

#define randomize() srandom((int)time(NULL))
/* initialize random number generator */

#define setSeed(N) srandom(N)

#define random01() ((double) random() / 2147483648.)
/* gives a floating point random number in [0,1) */

#define rnd(num) ((int) (random01() * (num)))
/* gives an integer random number between 0 and num - 1, inclusive */

// classes used in many places

class OptionLine;
class BinCoefTable;
class Move;
class Neighborhood;
class DataProvider;
class DistanceDataProvider;
class CoveringDataProvider;
class PatchMatrix;
class CoverMatrix;
class DistanceMatrix;
class CoverPenalty;
class DistancePenalty;
class Solution;
class CoverSolution;
class DistanceSolution;
class Stopwatch;
class Attribute;
class TabuList;
class Warning;
class Tenure;

// global variables

extern OptionLine *options;
extern int verbose;
extern Warning *warning;

#endif
