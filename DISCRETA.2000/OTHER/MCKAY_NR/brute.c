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
FILE *sol_txt;

static int lambda = 0;
static equation eqn[MAXEQN];
int **lo, **hi;
static int lorhs[MAXEQN], hirhs[MAXEQN], neqn[MAXEQN];
int nb_sol = 0;
long stopafter;
long stoploops;

int main(int argc, char **argv)
{
	int numvar,numeqn;
	int i,j,k;
	term *e;
	char zeile[100000], *fname;

	stoploops = 0;
	stopafter = 0;
	do {
		gets(zeile);  
      if (strstr(zeile,"% stopafter")!=NULL) {
			sscanf(zeile,"%% stopafter %ld",&stopafter);
		}
      if (strstr(zeile,"% stoploops")!=NULL) {
			sscanf(zeile,"%% stoploops %ld",&stoploops);
		}
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

printf("New version \n");
{
	int ii;
	
	lo=(int**)calloc(numvar,sizeof(int*));
	for (i=0;i<numvar;i++) lo[i]=(int*)calloc(numvar,sizeof(int));
	hi=(int**)calloc(numvar,sizeof(int*));
	for (i=0;i<numvar;i++) hi[i]=(int*)calloc(numvar,sizeof(int));

	for (i = 0; i < numeqn; i++) {
		e = (term*) eqn[i]; 
		k = 0;
	    
		for (j = 0; j < numvar; ++j) {
			scanf("%d", &ii);
			if (ii != 0) {
				e[k].var = j;
				e[k].coeff = ii;
				++k;
			}
		}
		neqn[i] = k;
/*
   \|lorhs| and  \|hirhs| are both equal to $\lambda$, because
   we want to solve equations, not inequations.
*/   
		lorhs[i] = lambda;
		hirhs[i] = lambda;
	} 
}

/* 
   We are searching for $\{0,1\}$-vectors. So the
   upper and lower bounds on the variables are $1$ respectively $0$.
*/
/*	for (i = 0; i< numvar; i++) { */
		for (j = 0; j < numvar; ++j)	{
			lo[0][j] = 0;
			hi[0][j] = 1;
		}
/*	} */

/* preprocessing phase:
   \|varprune| tries to find variables which are already fixed.
   \|divideeqn| reduces the equations, so that the gcd of all
   coefficients of each equation equals $1$.
*/
/* 	fprintf(out_txt,"orginal:\n"); 
	puteqns(lo[0],hi[0],numvar,lorhs,hirhs,eqn,neqn,numeqn);

 */	
 	varprune(lo[0],hi[0],lorhs,hirhs,eqn,neqn,numeqn);

/* 	fprintf(out_txt,"after varprune:\n");
	puteqns(lo[0],hi[0],numvar,lorhs,hirhs,eqn,neqn,numeqn);
 */
	if ((i = divideeqns(lorhs,hirhs,eqn,neqn,numeqn)) >= 0)	{
	    printf("division problem, equation %d\n",i);
	    exit(0);
	}
	fflush(stdout);

/* Now the solving routine starts. */

	sol_txt = fopen("solutions", "w");
	printf("Start of possolve numvar %d  numeqn %d\n",numvar,numeqn);
	possolve(lo,hi,eqn,lorhs,hirhs,neqn,numeqn,numvar);
	printf("Finito :-) Number of solutions: %ld\n\n", nb_sol);
	fflush(stdout);
	fclose(out_txt);
	fclose(sol_txt);
	exit(0); 
}

/* 
\subsection{solproc}
writes out solution vector.
*/
void solproc(int *x, int n)
{
	int i;

	for (i = 0; i < n; ++i) printf("%d",x[i]);
	printf(" l = %ld\n", lambda);
	fflush(stdout);

	for (i = 0; i < n; ++i) fprintf(sol_txt,"%d",x[i]);
	fprintf(sol_txt,"\n");
	fflush(sol_txt);

	nb_sol++;
	if ((stopafter>0)&&(nb_sol>=stopafter)) {
		fclose(sol_txt);
		exit(0);
	}
}


/* \endc */
