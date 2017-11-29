/* 
\input cnoweb
\title{general purpose linear diophantine (in)equation solver}
\job{brute.c}
\synopsis{brute.c : solve Kramer-Mesner equation as in Laue's paper.

The program was written by Brendan McKay.
It was ported to ANSI C by Alfred Wassermann, 1. March 1995.

It solves linear diophantine (in)equations with nonnegative coefficients
by a recursive brute force method.

It is possible to fix variables in the input file.}

\section{The main program}
Input format:

\item{(1)} Number of columns of the Kramer-Mesner-matrix

\item{(2)} Number of rows of the Kramer-Mesner-matrix

\item{(3)} $\lambda$

\item{(4)} The Kramer-Mesner-matrix without spaces.

\item{(5)} eventually fixed variables. First the number of fixed variables,
      then the pairs (number of -- variable value) e. g. 5 1.
This means the variable number 5 has the value 1.

*/
#include "brute.h"

/* 
\subsection{Main program} 

   There are some arrays which are declared \"globally".
   Maybe it's not necessary, but it works.
*/
  
FILE *out_txt;

static int lambda = 0;
static equation eqn[MAXEQN];
static int lo[MAXVAR],hi[MAXVAR];
static int lorhs[MAXEQN],hirhs[MAXEQN],neqn[MAXEQN];
int nb_sol = 0;

int main(int argc, char **argv)
{
	int numvar,numeqn;
	int i,j,k,nv,val,numstarts;
	char s[2];
	term *e;
	char zeile[80], *fname;

	do {
		gets(zeile);  
	}
	while (zeile[0]=='%');
	
/* input of the number of variables, number of equations and $\lambda$.
*/
	sscanf(zeile,"%d%d%d", &numeqn, &numvar);

	fname = argv[1];
	if (argc < 2) {
		printf("wrong number of args !\n");
		fflush(stdout);
		return 0;
		}
	lambda = atoi(argv[2]);
	printf("lambda = %d\n", lambda);
	fflush(stdout);
	out_txt = fopen(fname, "w");
	
#if 0
	if (scanf("%d%d%d", &numvar, &numeqn, &lambda) != 3) {
               fprintf(out_txt,">E Wrong input file format \n");
               exit(1);
	}
#endif

	if (numvar > MAXVAR || numeqn > MAXEQN) {
		fprintf(out_txt,">E increase MAXVAR or MAXEQN (two files)\n");
		printf(">E increase MAXVAR or MAXEQN (two files)\n");
		fflush(stdout);
		exit(1);
		}

/* input of the matrix on the left hand side.
   The right hand side vector consists only of $\lambda$ in
   every coordinate.
   
   Only coefficients different from zero are stored.
*/

	{
	int ii;
	
	for (i = 0; i < numeqn; i++) {
		e = (term*) eqn[i]; 
		k = 0;
	    
		for (j = 0; j < numvar; ++j) {

#if 0
	    if (scanf("%1s",&s[0]) != 1 || s[0] < '0' || s[0] > '9') { 
	    	    fprintf(out_txt,">E error in equation %d\n",i);
		    exit(1);
	        }
#endif
			scanf("%d", &ii);
	
			if (ii /* s[0] */ != 0 /* '0' */) {
				e[k].var = j;
				e[k].coeff = ii /*  (int) (s[0] - '0') */;
				++k;
				}
			} /* next j */
		neqn[i] = k;
/*
   \|lorhs| and  \|hirhs| are both equal to $\lambda$, because
   we want to solve equations, not inequations.
*/   
		lorhs[i] = lambda;
		hirhs[i] = lambda;
		} /* next i*/

	}

/* 
   We are searching for $\{0,1\}$-vectors. So the
   upper and lower bounds on the variables are $1$ respectively $0$.
*/
	for (j = 0; j < numvar; ++j)
	{
	    lo[j] = 0;
	    hi[j] = 1;
	}

#if 0
/* Additional input. Here we can fix certain variables. */
	if (scanf("%d",&numstarts) == 1)
	{
	    for (j = 0; j < numstarts; ++j) {
		if (scanf("%d%d",&nv,&val) != 2) {
		    fprintf(out_txt,">E can't read starter %d\n",j);
		    exit(1);
		}
	  	if (nv < 0 || nv >= numvar) {
		    fprintf(out_txt,">E starter variable out of range\n");
		    exit(1);
		}
		lo[nv] = val;
		hi[nv] = val;
	    }
	}
#endif
	numstarts = 0;
	
/* preprocessing phase:
   \|varprune| tries to find variables which are already fixed.
   \|divideeqn| reduces the equations, so that the gcd of all
   coefficients of each equation equals $1$.
*/
	fprintf(out_txt,"orginal:\n"); 
	puteqns(lo,hi,numvar,lorhs,hirhs,eqn,neqn,numeqn);

	varprune(lo,hi,lorhs,hirhs,eqn,neqn,numeqn);

	fprintf(out_txt,"after varprune:\n");
	puteqns(lo,hi,numvar,lorhs,hirhs,eqn,neqn,numeqn);

	if ((i = divideeqns(lorhs,hirhs,eqn,neqn,numeqn)) >= 0)
	{
	    printf("division problem, equation %d\n",i);
	    exit(0);
	}
	fflush(stdout);

/* Now the solving routine starts. */
	possolve(lo,hi,eqn,lorhs,hirhs,neqn,numeqn,numvar);
	printf("Finito :-) %ld\n\n", nb_sol);
	fflush(stdout);
	fclose(out_txt);
        return 1; 
}

/* 
\subsection{solproc}
writes out solution vector.
*/
void solproc(int *x, int n)
{
	int i;

	fflush(out_txt);
	/* printf("S"); */
	for (i = 0; i < n; ++i)
	    printf("%d ",x[i]);
	printf("l = %ld\n", lambda);
	fflush(stdout);
	nb_sol++;
}


/* \endc */
