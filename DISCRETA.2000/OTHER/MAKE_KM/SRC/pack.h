/*   pack.h
**
**   This file contains definitions common to all modules (or many of
**   them). This file should be included in *every* module.
**
**
**  Brute force Kramer-Mesner matrix generation.
**
**  (c) 1997,1999 Kari J. Nurmela
*/

#ifndef _pack_h_
#define _pack_h_

#include <stdio.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/time.h>
#include <math.h>
#include <values.h>

/***** type definitions *****/

typedef unsigned cword;
#define CWORD_BITS (8*sizeof(cword))
/* codewords are of this type */
/* cword MUST be unsigned */

typedef unsigned short orbitSize_t;

#define bitCount(X,C) {tmpcword = X; C = bits8[tmpcword % 256]; \
tmpcword /= 256;\
C += bits8[tmpcword % 256]; tmpcword /= 256; C += bits8[tmpcword % 256]; \
tmpcword /= 256;\
C += bits8[tmpcword];}
extern cword tmpcword;

/***** static table sizes *****/

#define MAX_ORBIT_LEN 2000000

/***** global variables *****/

extern cword *ctable;
extern int ctableSize;
extern cword orbit[];

/***** functions *****/

extern void printGenerators(FILE *, int);
extern int generateOrbitAndSetCtableBits(cword, cword*, int);
extern int generateOrbit(cword, cword*, int);

#endif
