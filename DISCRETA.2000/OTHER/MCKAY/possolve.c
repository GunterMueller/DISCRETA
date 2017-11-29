/*
\input cnoweb
\title{Linear Diophantine Equation Solver}
\job{possolve.c}
\synopsis{\|possolve| finds solutions for linear
diophantine (in)equations with nonnegative coefficients 
by a backtracking method.

It was written by Brendan McKay and changed by Alfred Wassermann.
It is especially tuned for solving $\{0,1\}$-problems.}


Debugging output:

If \|verbose| $ = 0$: No debugging output. 

If \|verbose| $ = 1$: No debugging output.

If \|verbose| $ = 2$: Branch and split information.

If \|verbose| $ = 3$: detailed information about equations and variables.

\section{The header file \|brute.h|}
\includefile{brute.h}

\section{possolve}
\|possolve| finds  solutions of linear diophantine (in)equations
 with nonnegative coefficients.
 
 First the bounds on the variables are restricted as far as possible,
 then the recursive procedure \|solve| is started.
 
*/
#include "brute.h"

/* There are two global arrays --- \|unitcoeffs| and \|active| ---
   which are declared in all subroutines of the file \|possolve.c|.
   
   \|unitcoeffs[i]| is true if the $i$th equation is a $\{0,1\}$-equation. 
   \|active| denotes all equations which are not pruned.
*/
static boolean unitcoeffs[MAXEQN];
static boolean active[MAXEQN];


/* \subsection{possolve} 

  solve (in)equations with positive integer coefficients.

  numvar = number of variables
  numeqn = number of equations
  
     (The bounds MAXVAR and MAXEQN are defined above.  MAXVAR must be
      the same in any program calling this procedure.)

  \|lo[i]|,\|hi[i]| $=$ initial bounds on variable $i$.
  (They should be $\geq 0$.)

  \|aeqn[j]| = the $j$-th (in)equation.  
  (It is a list of (variable,coefficient) where all the coefficients
      must be strictly positive.  The number of terms is \|aneqn[j]|.)
   
  \|alorsh[j]|, \|ahirhs[j]| = bounds on the value of the $j$-th (in)equation.

     (So, the statement is \|alorhs[j]| $\leq$ \|aeqn[j]| $\leq$ \|ahirhs[j]|.
      It must be that these bounds are $\geq 0$.)

  \|asolproc(x,numvar) int x[],numvar;| = a procedure to call for each solution.
  (The solution is \|x[0],x[1],...,x[numvar-1]|.)

  Comments: The code is tuned for short equations, and likes equations
     with 0-1 coefficients especially.

  Brendan McKay, Feb 1995.
*/

int rekurs = 0;

void possolve(int *lo, int *hi, equation *eqn, int *lorhs, int *hirhs, int *neqn,
int numeqn, int numvar)
{
        register int i,j;
        register term *e;
	boolean hopeless;

        if (numeqn > MAXEQN || numvar > MAXVAR) {
            fprintf(out_txt,">E intsolve: limits exceeded\n");
            exit(1);
        }
/* First step: 
   If for one equation there is no coefficient different from
   zero, this equation is \"not" active. 
   Further mark in \|unitcoeffs| the equations with coefficients
   solely $\in\{0,1\}$.
*/
	hopeless = FALSE;
        for (i = 0; i < numeqn; ++i) {
	    if (neqn[i] == 0) {
		active[i] = FALSE;
		if (lorhs[i] > 0 || hirhs[i] < 0) hopeless = TRUE;
	    }
	    else
		active[i] = TRUE;
            e = (term*) eqn[i];
            for (j = neqn[i]; --j >= 0;) if (e[j].coeff != 1) break;
            unitcoeffs[i] = (j < 0);
        }

/* More output. */
#if VERBOSE > 2
        	puteqns(lo,hi,numvar,lorhs,hirhs,eqn,neqn,numeqn);
#endif

/* Finally the recursion starts. */
        if (!hopeless)
	        solve(0,lo,hi,active,numvar,lorhs,hirhs,eqn,neqn,numeqn);
	return;
}

/* \subsection{subtract} 
  subtract equation e1 from e2 if possible  ---
 return success.
 It is used in function \|pruneqn|.
*/
boolean subtract(term *e1, int l1, int lors1, int hirs1, term *e2, int *pl2, int *plors2, int *phirs2)
{                              
        register int i,j,k;
        term e1i;
        int l2,factor,minfactor;

/* First test if subtraction is possible. */
        minfactor = 999999;
        l2 = *pl2;
        if (l1 > l2 || hirs1 > lors1) return FALSE;

/* Second test if subtraction is possible. */
        j = 0;
        for (i = 0; i < l1; ++i) {
            e1i = e1[i];
            for (; j < l2 && e2[j].var != e1i.var; ++j) {}
            if (j == l2 || e2[j].coeff < e1i.coeff) return FALSE;
            factor = e2[j].coeff / e1i.coeff;
            if (factor < minfactor) minfactor = factor;
        }
        
/* Do subtraction */
        k = 0;
        for (i = j = 0; i < l1; ++i, ++j)
        {
            e1i = e1[i];
            for (; j < l2 && e2[j].var != e1i.var; ++j) e2[k++] = e2[j];
            if (j < l2 && e2[j].coeff > minfactor*e1i.coeff) {
                e2[k].var = e2[j].var;
                e2[k].coeff = e2[j].coeff - minfactor*e1i.coeff;
                ++k;
            }
        }
        for (; j < l2; ++j) e2[k++] = e2[j];

        *pl2 = k;
        *plors2 -= minfactor*lors1;
        *phirs2 -= minfactor*hirs1;

/*  end of subtraction. */
        return TRUE;
}

/* \subsection{pruneqn --- remove equations} 
 prune equations by subtraction.
*/
void pruneqn(int *lorhs, int *hirhs, equation *eqn, int *neqn, int numeqn)    
{
        boolean ok,done[MAXEQN];
        register int i,j;

/* \|done| is always \|FALSE|. Why? */

        for (i = 0; i < numeqn; ++i) done[i] = FALSE;

        do {
            ok = TRUE;
            for (i = 0; i < numeqn; ++i)
            if (!done[i] && neqn[i] > 0) {
                for (j = 0; j < numeqn; ++j)
                if (i != j && subtract(eqn[i],neqn[i],lorhs[i],hirhs[i],
                                       eqn[j],&neqn[j],&lorhs[j],&hirhs[j])) {
                    ok = FALSE;
                    done[j] = FALSE;
                }
            }
        } while (!ok);
        return;
}

/* 
  \subsection{varprune --- remove variables} 
   Try to remove free variables by testing if variables are already
   fixed. This is the case if the lower bound on a variable is equal to its
   upper bound.
*/

void varprune(int *lo,int *hi,int *lorhs,int *hirhs,equation *eqn,int *neqn,int numeqn)
{
        register int i,j,sum,len;
        register term *e;

        for (j = 0; j < numeqn; ++j) {
            e = (term*)eqn[j];
            len = neqn[j];
            sum = 0;
/* simple test whether the lower bound of a variable meets its upper bound. */
            for (i = 0; i < len;) {
                if (lo[e[i].var] == hi[e[i].var]) {
                    sum += e[i].coeff*lo[e[i].var];
                    e[i] = e[--len];
                }
                else
                    ++i;
            }
            lorhs[j] -= sum;
            hirhs[j] -= sum;
            neqn[j] = len;
        }
}

/* This function seems to be an relict from an older version. */

void nul() {}
/* 
\subsection{Output routine --- puteqns} 
write equations in dildo style to out_txt. 
*/
void puteqns(int *lo,int *hi,int numvar,int *lorhs,int *hirhs,equation *eqn,int *neqn,int numeqn)
{
        int i,j;
        term *e;

/*  First the lower and upper bounds on the variable are written. */
        fprintf(out_txt,"CON\n");
        for (i = 0; i < numvar; ++i) {
            if (i > 0) fprintf(out_txt,",\n");
            if (lo[i] > 0) fprintf(out_txt,"x%d >= %d, ",i,lo[i]);
            fprintf(out_txt,"x%d <= %d",i,hi[i]);
        }
/* Now the equations are written. */
        for (i = 0; i < numeqn; ++i) {
            e = (term*) eqn[i];
            fprintf(out_txt,",\n");
            for (j = 0; j < neqn[i]; ++j) {
                if (j % 16 == 15) fprintf(out_txt,"\n  ");
                fprintf(out_txt,"%s%d*x%d",(j == 0 ? "" : "+"),e[j].coeff,e[j].var);
            }
            fprintf(out_txt,"%s%d",lorhs[i] < hirhs[i] ? " >= " : " = ",lorhs[i]);

            if (lorhs[i] == hirhs[i]) continue;

/* Here is the output of inequations. */
            fprintf(out_txt,",\n");
            for (j = 0; j < neqn[i]; ++j) {
                if (j % 16 == 15) fprintf(out_txt,"\n  ");
                fprintf(out_txt,"%s%d*x%d",(j == 0 ? "" : "+"),
                               e[j].coeff,e[j].var);
            }
            fprintf(out_txt,"%s%d"," <= ",hirhs[i]);
        }
        fprintf(out_txt,";\n");
}

/* 
 \subsection{divideeqns} 
 take out common factors, return bad eqn number.
 It is only used in the main program.
*/
int divideeqns(int *lorhs, int *hirhs, equation *eqn, int *neqn, int numeqn)
{
	register int i,j,g,len;
	register term *e;

	for (j = 0; j < numeqn; ++j)
        {
            e = (term*) eqn[j];
            len = neqn[j];
	    if (len == 0) continue;

	    g = e[0].coeff;
	    i = 1;
	    for (i = 1; i < len && g > 1; ++i) g = gcd(g,e[i].coeff);
/* $g = \gcd$ of all coefficients of the left hand side of equation $i$. 
   If $g=1$ step to the next equation.
*/
	    if (g == 1) continue;
	    for (i = 0; i < len; ++i) e[i].coeff /= g;

	    lorhs[j] = lorhs[j] < 0 ? 0 : (lorhs[j] + g - 1) / g;
	    hirhs[j] = hirhs[j] < 0 ? -1 : hirhs[j] / g;

/*  Write some information.*/
	    fprintf(out_txt, "eqn %d: g=%d lorhs=%d hirhs=%d\n",j,g,lorhs[j],hirhs[j]);

	    if (lorhs[j] > hirhs[j]) return j;
	}

	return -1;
}

/* 
\subsection{gcd}
used in \|divideeqns|.
*/ 
int gcd(int n1,int n2)
{
        register int a,b,c;
 
        if (n1 > n2) {
            a = n1; b = n2;
        }
        else {
            a = n2; b = n1;
        }
 
        while ((c = a % b) > 0) {
            a = b; b = c;
        }
 
        return(b);
}

/*
\section{solve --- brute force recursion}
This procedure is called recursively.
*/
void solve(int level,int *alo,int *ahi,boolean *aactive,int numvar,int *lorhs,int *hirhs,equation *eqn,int *neqn,int numeqn)
{
        register int i,j,losum,hisum,eic,eiv,lx,hx;
        int nfree,lo[MAXVAR],hi[MAXVAR],ok;
        int xlo,xhi,isplit,mindiff;
        term *e;
        boolean active[MAXEQN];

/* Some debugging and informational output. */
#if VERBOSE > 1
        	int nacc;
#if VERBOSE > 3
			fprintf(out_txt," SOLVE level %d\n",level);
#endif
        	if (ticks++ % INTERVAL == 0) {
	            nacc = 0;
        	 /* *< 
	                for (i = 0; i < level && i < 15; ++i)
        	             fprintf(out_txt," %d/%d",range[i],split[i]);
	                fprintf(out_txt,"\n"); 
        	 >*  */
	            fprintf(out_txt,"Number of active equations: ");
        	    for (i = 0; i < level; ++i) fprintf(out_txt,"%d",branch[i]); 
	            for (i = numeqn; --i >= 0;) if (aactive[i]) ++nacc; 
        	    fprintf(out_txt,":%d\n",nacc);
	        }
#endif

/*  \|lo|, \|hi| and \|active| are local arrays. */

        for (i = 0; i < numvar; ++i) {
            lo[i] = alo[i];
            hi[i] = ahi[i];
        }
        for (i = 0; i < numeqn; ++i) active[i] = aactive[i];

/* The following line seems to be a relict from another proplem. */
/* *<        if (level == 8 && lo[22] == 1 && lo[19] == 1) nul(); >* */

/* The following loop through the equations tries to restrict the lower
   and upper bounds on the variables.
   We only have to handle active equations.
*/
        ok = 0;
	/* ok = number of equations that have been 
	 * checked since the last change of any 
	 * lo[] or hi[] was made. 
	 * the aim is to check all equations once 
	 * without reducing hi[] - lo[]; 
	 * so that we have reached stability */
        for (j = 0; ok < numeqn; j = (j == numeqn-1 ? 0 : j+1)) {
		/* j is the next equation to check;
		 * j wraps around all equation indices */
            ++ok;
            if (active[j]) {
			/* we check equation j: */

/* We distinguish two cases: First if there are only $\{0,1\}$-coefficients. */
                if (unitcoeffs[j]) {
			e = (term*)eqn[j];
/* The lower and upper bounds on the variables (belonging to
   nonzero coefficients) are summed up in \|losum| and \|hisum|.
*/   
			losum = 0;
			hisum = 0;
			for (i = neqn[j]; --i >= 0;) {
				losum += lo[e[i].var];
					/* lowest possible rhs */
				hisum += hi[e[i].var];
					/* highest possible rhs */
				}
			if (losum > hirhs[j] || hisum < lorhs[j])
				return;

/* 
   If possible the lower or upper bounds on the variables are further
   restricted.
	count the number of remaining free variables in nfree. 
*/
			nfree = 0;
			for (i = neqn[j]; --i >= 0;) {
				eiv = e[i].var;
				hx = hi[eiv];
				lx = lo[eiv];
				if (hx != lx) {
					xlo = lorhs[j] + hx - hisum;
					xhi = hirhs[j] + lx - losum;
					if (xlo > lx) {
						lo[eiv] = xlo;
						ok = 0;
						/* a change was made;
						 * loop through all 
						 * equations again. */
						}
					if (xhi < hx) {
						hi[eiv] = xhi;
						ok = 0;
						}
					if (lo[eiv] != hi[eiv])
						++nfree;
					} /* if */
				} /* next i */
			} /* if (unitcoeffs[j]) */
/* Now the slightly more complicated case if there are coefficents greater
   than $1$.
   If the lower bound of a variable becomes greater than its upper bound,
   the procedure is stopped at once.
*/
		else {
			e = (term*) eqn[j];
/* Again the lower and upper bounds on the variables (belonging to
   nonzero coefficients) are summed up in \|losum| and \|hisum|. 
*/   
			losum = 0;
			hisum = 0;
			for (i = neqn[j]; --i >= 0;) {
				losum += e[i].coeff * lo[e[i].var];
				hisum += e[i].coeff * hi[e[i].var];
				}
			if (losum > hirhs[j] || hisum < lorhs[j])
				return;

/* 
   And if possible the lower or upper bounds on the variables are further
   restricted.
*/
			nfree = 0;
			for (i = neqn[j]; --i >= 0;) {
				if (hi[e[i].var] != lo[e[i].var]) {
					eic = e[i].coeff;
					eiv = e[i].var;
					hx = eic * hi[eiv];
					lx = eic * lo[eiv];
					xlo = lorhs[j] + hx - hisum;
					xhi = hirhs[j] + lx - losum;
					if (xlo > lx) {
						lo[eiv] = (xlo + eic - 1) / eic;
						if (lo[eiv] > hi[eiv])
							return;
						ok = 0;
						}
					if (xhi < hx) {
						hi[eiv] = xhi / eic;
						if (lo[eiv] > hi[eiv])
							return;
						ok = 0;
						}
					if (lo[eiv] != hi[eiv])
						++nfree;
					} /* if */
				} /* next i */
			} /* else */
/* Now hopefully the variables are in each case further restricted.
   The equation becomes inactive if

   \item{(1)} there are no free variables in this equation \"and"
   
   \item{(2)} if it was not possible to further restrict the variables in the
         last try.
*/         
                if (ok > 0 && nfree == 0) active[j] = FALSE;
            }
        }

/* Here comes the searching part. \|mindiff| gives the smallest
   difference between lower and upper bound of a variable. \|isplit|
   is the corresponding variable number.
*/
        mindiff = MAXB+1;
        isplit = -1;
        for (i = 0; i < numvar; ++i)
        if (hi[i] != lo[i] && hi[i] - lo[i] < mindiff) {
            isplit = i;
            mindiff = hi[i] - lo[i];
        }

/* If \|isplit| $< 0$ we have found a solution. Otherwise we try to
  delete variables by \|varprune| and we try to delete equations
  by \|pruneqn|.
*/
        if (isplit < 0)  (solproc)(lo,numvar);
        else {
            if (level == 0) {
                varprune(lo,hi,lorhs,hirhs,eqn,neqn,numeqn);
                
#if VERBOSE > 2
	                puteqns(lo,hi,numvar,lorhs,hirhs,eqn,neqn,numeqn);
#endif

		
                pruneqn(lorhs,hirhs,eqn,neqn,numeqn);

#if VERBOSE > 2
	                puteqns(lo,hi,numvar,lorhs,hirhs,eqn,neqn,numeqn);
#endif
            }

/* 
   Finally we start the recursion.
   the variable with the number \|isplit| runs from
   \|lo[isplit]| to \|hi[isplit]| and for each step we call \|solve|
   with this fixed value.
 
   \|branch|, \|split| and \|range| are collected for debugging purposes. 
*/
            for (i = 0; i <= mindiff; ++i)
            {
                hi[isplit] = lo[isplit];
#if VERBOSE > 1
                	branch[level] = lo[isplit];
	                split[level] = isplit;
        	        range[level] = mindiff+1;
#endif
                solve(level+1,lo,hi,active,numvar,lorhs,hirhs,eqn,neqn,numeqn);
                ++lo[isplit];
            }
        }
}

/* \endc */
