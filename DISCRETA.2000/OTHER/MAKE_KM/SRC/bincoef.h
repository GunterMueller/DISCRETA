/*  Brute force Kramer-Mesner matrix generation.
**
**  (c) 1997,1999 Kari J. Nurmela
**
**   bincoef.h
**
**   This file contains definitions needed when calculating the binomial
**   coefficient table and when using that table. 
**
*/


#ifndef _bincoef_h_
#define _bincoef_h_

/***** type definitions *****/

#define maxv 40
/* binomial coefficients are tried to calculate up to binCoef[maxv][?]
 * (overflow is checked and the program is not halted)
 */

typedef unsigned binCoefType;
/* this type should be able to contain binomial coefficients needed */

/***** macros *****/

#define bincoef(A,B) (binCoef[A][B])
#define overflowBinCoef(v,k) (binCoef[(v)][(k)] == 0)

/***** global variables *****/

extern binCoefType binCoef[maxv + 1][maxv + 2];

/***** functions *****/

void calculateBinCoefs(void);
extern float bincoefF(int, int);

#endif
