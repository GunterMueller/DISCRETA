/* % $Header: /usr/local/cvsroot/designs/llll.h,v 1.9 1997/12/03 17:40:26 alfred Exp $ */
#ifndef _LLLL_H
#define _LLLL_H

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include "rueck_werte.h"

/* #define GIVENS 1 */
#define DOUBLE double
#define SQRT sqrt
#define VERBOSE 1
#define EPSILON 0.00001
#define LLLCONST_LOW 0.75
#define LLLCONST_HIGH 0.99

struct coe {
	long c;
	int p;
};
#define COEFF struct coe

extern int s_sel;
extern int inequalities;
extern int slackvar;
extern long stopafter;
extern long stoploops;

extern long maximal_number_of_blocks;
extern long midpoint;
extern long grouporder;
extern long lambda;

extern int lll (COEFF **b, int s, int z, int KMz, long v1, long v2, long c1, long nom, long denom, int flag);
extern DOUBLE scalarproductfp (DOUBLE *v, DOUBLE *w, int n);
extern long scalarproductl (COEFF *v, COEFF *w/*, int n*/);
extern DOUBLE scalarproductlfp (COEFF *v, COEFF *w/*, int n*/);
extern long norml (COEFF *v/*, int n*/);
extern DOUBLE normfp (COEFF *v/*, int n*/);
extern int lllfp (COEFF **b, DOUBLE **mu, DOUBLE *c, DOUBLE *N, /*DOUBLE **bs,*/ int start, int s, int z, DOUBLE delta, int KMz, long v1, long v2, long c1, long nom, long denom, int flag);
extern int bkz (COEFF **b, int s, int z, DOUBLE delta, int beta, int KMz, long v1, long v2, long c1, long nom, long denom, int flag, int p );
extern DOUBLE round (DOUBLE r);
extern void ausgabe (COEFF **b, int s, int z);
extern void vectorausgabe (COEFF *v, int z);
extern double analyse (COEFF **b, int s, int z);
extern DOUBLE enumerate (DOUBLE **mu, DOUBLE *c, long *u, int s, int j, int k, int p);
extern DOUBLE smallenumerate (DOUBLE **mu, DOUBLE *c, long *u, int s, int z, int j, int k);
extern int lllalloc(DOUBLE ***mu, DOUBLE **c, DOUBLE **N, /*DOUBLE ***bs,*/ int s, int z);
extern int lllfree(DOUBLE **mu, DOUBLE *c, DOUBLE *N, /*DOUBLE **bs,*/ int s);
extern void putlattice(COEFF **b, int s, int z, int KMz, int KMs);
extern void getlattice(COEFF **b, int *s, int *z, int flag);
extern void cutlattice(COEFF **b, int *s, int *z, int KMz);
extern void coeffinit(COEFF *v, int z);

extern int solutiontest(COEFF *b, int z, int KMz, long v1, long v2, long c1, long nom, long denom, int flag);
extern DOUBLE enuminf (COEFF **b, int s, int z, int KMs, int KMz, long v1, long v2, /*long c0, long c1,*/ long nom, long denom, int selection, int *permutation);

#endif /* _LLLL_H */
