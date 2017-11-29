#ifndef BRUTE_H
#define BRUTE_H

#define EXTDEFS
/* *< #include "des.h" >* */
#include <stdio.h>
#include <math.h>

/* bigger gets more diagnostic output */
#define VERBOSE 0  

/* Some constants needed in brute.c and possolve.c. */
#define MAXB   1000 /* 1000 */
#define MAXEQN 100 /* 100 */
#define MAXVAR 250 /* 127 */

#define INT short
#define FALSE 0
#define TRUE  1
#define boolean int

/* The main structure: \|term| which contains the number of th coefficient
   and its value. 
   \|equation| is an array of \|term|s.
*/
typedef struct {int var,coeff;} term;
typedef term equation[MAXVAR+1];

#if VERBOSE > 1
static int range[MAXEQN+1],split[MAXEQN+1],branch[MAXEQN+1];
static int ticks = 0;
#define INTERVAL 10000
#endif

extern FILE *out_txt;

/* Here are the function prototypes. */
extern boolean subtract(term *e1, int l1, int lors1, int hirs1, term *e2, int *pl2, int *plors2, int *phirs2);
extern void pruneqn(int *lorhs, int *hirhs, equation *eqn, int *neqn, int numeqn);   
extern void varprune(int *lo, int *hi,int *lorhs, int *hirhs, equation *eqn, int *neqn, int numeqn);
extern void nul();
extern void solve(int level, int *alo, int *ahi, boolean *aactive, int numvar,
int *lorhs, int *hirhs, equation *eqn, int *neqn, int numeqn);
extern void solproc(int *lo, int numvar);
extern void puteqns(int *lo,int *hi,int numvar,int *lorhs,int *hirhs,equation *eqn,int *neqn,int numeqn);
extern void possolve(int *lo, int *hi, equation *eqn, int *lorhs, int *hirhs, int *neqn, int numeqn, int numvar);

extern int divideeqns(int *lorhs, int *hirhs, equation *eqn, int *neqn, int numeqn);
extern int gcd(int n1,int n2);

#endif
