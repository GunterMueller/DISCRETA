/* 
% $Header: /usr/local/cvsroot/designs/enumall.c,v 1.10 1997/12/03 17:40:24 alfred Exp $
\section{Explicit enumerate of (all) solutions}

Changed: 30.9.1996,
	29.10.1996.
	30.10.1996: Additional compilerflag: GIVENS (set in \|llll.h|)

Original from Kaib/Ritter.
*/

#include "llll.h"
long nosolutions;

#if defined(PRINTSOLUTION)
	static FILE *fp;

void open_solution_file() {
	fp = fopen("solutions", "w");
}

void close_solution_file() {
	fclose(fp);
}
#endif

/* These functions are separated from the main code to measure
   their performance. If the compiler uses inline functions
   it does not slow down the computation.
*/   
DOUBLE compute_y(DOUBLE **mu, long *us, int t, int tmax) {
	int i;
	DOUBLE dum;

	for (i=t+1,dum=0.0;i<=tmax;i++) dum += mu[i][t]*us[i];  
	return dum;
}

void compute_w(DOUBLE **w, DOUBLE **bd, DOUBLE dum, int t, int z) {
	int l;
	
	for (l=z-1;l>=0;l--) {
		w[t][l] = w[t+1][l] + dum*bd[t][l]; 
	}		
	return;
}

/* Some prototypes */
DOUBLE loggamma(DOUBLE x);
int exacttest(DOUBLE *v, int z, DOUBLE Fq);
int prune0(DOUBLE li, DOUBLE re);
int prune(DOUBLE *w, DOUBLE cs, int z, DOUBLE Fq);
int pruneprimal(DOUBLE *w, int t, int z, DOUBLE Fq, int *firstlist, int *firstp);
int onlinesolutiontest(DOUBLE *w, int z, long c1, long nom, long denom, int selection, int *permutation);
void gramschmidt(COEFF **b, int k, int z, DOUBLE **mu, DOUBLE **bd, DOUBLE *c, DOUBLE *N, long solutionnorm);
void givens(COEFF **b, int k, int z, DOUBLE **mu, DOUBLE **bd, DOUBLE *c, DOUBLE *N, long solutionnorm);

/* \section{Enumeration}

The main enumeration routine. 
*/

DOUBLE enuminf (COEFF **b, int s, int z, int KMs, int KMz, long v1, long v2, /*long c0, long c1,*/ long nom, long denom, int selection, int *permutation)
{
	COEFF *swapvl;
	DOUBLE *y;
	DOUBLE *cs; 
	DOUBLE dum;
	long *us;
	int t,i,l,ii;
	int j,k;
	
	long *delta;
	long *d;
	long *eta;
	long *v;
	int *first, *firstlist, *firstp; 
	DOUBLE *N;
	DOUBLE **mu;
	DOUBLE *c;
	DOUBLE **w;
	DOUBLE **bd;
	DOUBLE Fd, Fq; 
	int tmax;
	long solutionnorm; 
	long loops = 0;

/*%--------------------------------------
*/
#if defined(PROBABILISTIC)
/* Blocks with blocksize lower than \|SCHNITT| are completely enumerated,
   otherwise the enumeration is pruned. */
	DOUBLE *gap;
	int p = 60;
	int SCHNITT = 10;
	DOUBLE help;
	
/* $2^{-p}$ is the probability that lattice points are lost. */
	DOUBLE pi = 3.141592653589793238462643383;
#endif	
/*%--------------------------------------
*/
	
	j = 0;
	k = s-1;

	lllalloc(&mu,&c,&N,/*&bd,*/s,z);

/* Output of the dimension of the solution space. */    
	printf("Dimension of solution space (k): %d compared to s-z+2: %d\n",k+1, KMs-KMz+2);
	fflush(stdout);             
	
/* Exit, if we don't have a basis of the whole kernel. */
	if (k+1 < KMs-KMz+2) {
		printf("The LLL algorithm didn't succeed in computing a basis of the kernel.\n");
		printf("Please increase c0 (the first parameter)!\n");
		return 0;
	}
	
/* allocate memory for \|us| */
     us=(long*)calloc(s+1,sizeof(long));
     if (us == NULL) return(zu_wenig_Speicherplatz);

/* allocate memory for \|cs| */
     cs=(DOUBLE*)calloc(s+1,sizeof(DOUBLE));
     if (cs == NULL) return(zu_wenig_Speicherplatz);

/* allocate memory for \|y| */
     y=(DOUBLE*)calloc(s+1,sizeof(DOUBLE));
     if (y == NULL) return(zu_wenig_Speicherplatz);

/* allocate memory for \|delta| */
     delta=(long*)calloc(s+1,sizeof(long));
     if (delta == NULL) return(zu_wenig_Speicherplatz);

/* allocate memory for \|d| */
     d=(long*)calloc(s+1,sizeof(long));
     if (d == NULL) return(zu_wenig_Speicherplatz);

/* allocate memory for \|first| */
     first=(int*)calloc(z,sizeof(int));
     if (first == NULL) return(zu_wenig_Speicherplatz);
     firstlist=(int*)calloc(s+z+1,sizeof(int));
     if (firstlist == NULL) return(zu_wenig_Speicherplatz);
     firstp=(int*)calloc(s+1,sizeof(int));
     if (firstp == NULL) return(zu_wenig_Speicherplatz);

/* allocate memory for \|eta| */
     eta=(long*)calloc(s+1,sizeof(long));
     if (eta == NULL) return(zu_wenig_Speicherplatz);

/* allocate memory for \|v| */
     v=(long*)calloc(s+1,sizeof(long));
     if (v == NULL) return(zu_wenig_Speicherplatz);

/* allocate memory for \|w| */
     w=(DOUBLE**)calloc(s+1,sizeof(DOUBLE*));
     if (w == NULL) return(zu_wenig_Speicherplatz);
     for (i=0;i<=s;i++) w[i]=(DOUBLE*)calloc(z,sizeof(DOUBLE));

/* allocate memory for \|bd| */
/*  Givens rotation */
#if defined(GIVENS)
     bd=(DOUBLE**)calloc(z,sizeof(DOUBLE*));
     if (bd == NULL) return(zu_wenig_Speicherplatz);
     for (i=0;i<z;i++) bd[i]=(DOUBLE*)calloc(z,sizeof(DOUBLE));
#else     
     bd=(DOUBLE**)calloc(s,sizeof(DOUBLE*));
     if (bd == NULL) return(zu_wenig_Speicherplatz);
     for (i=0;i<s;i++) bd[i]=(DOUBLE*)calloc(z,sizeof(DOUBLE));
#endif

/*%----------------------------------------------------------------*/
/* allocate memory for \|gap| */
#if defined(PROBABILISTIC)
     gap=(DOUBLE*)calloc(s+1,sizeof(DOUBLE));
     if (gap == NULL) return(zu_wenig_Speicherplatz);
#endif
/*%----------------------------------------------------------------*/


/* \subsection{Initiation}

of the arrays.
 */
 for (i=j;i<=k+1;i++) {
	cs[i] = 0.0;
	us[i] = 0;
	y[i] = 0.0;
	v[i] = 0;
	delta[i] = 0;
	eta[i] = 1;
	d[i] = 1;
	for (l=0;l<z;l++) w[i][l] = 0.0;
 }
	

/* \subsection{Reordering}

Reordering such that the columns with a nonzero entry in the last row
are at the end. 
*/

/* Original, before 22.10.96. */
/* *<       for (l=0;l<s/2;l++) {
	  if (labs(b[l][z].c)>EPSILON) {
	        swapvl = b[l];
		for (i=l+1;i<s;i++) b[i-1] = b[i];
		b[s-1] = swapvl;
		}
	}
>* */

/* 
	New version, since 22.10.96    
	Error occurs in rare cases: 27.6.1997.
	Replaced by simple sorting algorithm.
*/	
/* *<
	ii = s-1;
	for (l=s-2;l>=0;l--) {
		if (labs(b[l][z].c)>labs(b[ii][z].c)) {
	   	swapvl = b[l];
			for (i=l+1;i<=ii;i++) b[i-1] = b[i];
			b[ii] = swapvl;
			while ( (ii>0) && (labs(b[ii][z].c)>0)) ii--;
			if (ii==0) l= -1;
		}
	}
>* */
	 
/* simple sorting */

   for (ii=s-1;ii>0;ii--) {
	   for (l=ii-1;l>=0;l--) {
			if (labs(b[l][z].c)>labs(b[ii][z].c)) {
	     		swapvl = b[l];
				for (i=l+1;i<=ii;i++) b[i-1] = b[i];
				b[ii] = swapvl;
			}
		}
	}

/* \subsection{Setting of the upper bounds}

\|solutionnorm| is the maximum norm of a solution.
\|Fq| is the maximum norm of the solutions.
\|Fd| is the square of the $l_2$ norm of the solution vectors.

\|EPSILON| is a security buffer to avoid truncation due to rounding errors.
*/
    if (labs(v1)>labs(v2)) { solutionnorm = labs(v1); }
    else { solutionnorm = labs(v2); }

    Fq = (DOUBLE)solutionnorm; 
    Fd = (z*Fq*Fq)*(1.0+EPSILON); 

#if VERBOSE > 0 	
	printf("Fd: %f\n",(double)Fd);
#endif
/* Orthogonalization */
#if defined(GIVENS)
	givens(b,k,z,mu,bd,c,N,solutionnorm);  
#else	
	gramschmidt(b,k,z,mu,bd,c,N,solutionnorm);
#endif	

/* Since we only want to search in the affine subspace
with the last coefficent $\neq 0$ we have to test if we still are
searching in this affine subspace.
*/     
	for (l=0;l<z;l++) {
		for (i=0; i<=k; i++) if (labs(b[i][l+1].c) > 0) { 
			first[l] = i;
			break;
		}           
	}
	ii = 0;
	for (l=0;l<s;l++) {
		firstp[l] = ii;
		firstlist[ii] = 0;
		ii++;
		for (i=0;i<z;i++)
			if(first[i]==l) {
				firstlist[ii] = i;
				firstlist[firstp[l]]++;
				ii++;
			}
	}
	firstp[s] = ii;
	firstlist[ii] = 0; 

	nosolutions = 0;

	/* Setting the initial level. */
	t = first[z-1]-1;
	if (t<0) t = 0;
	tmax = t;
	us[t] = 1;
	v[t] = 1;

/*%---------------------------------------------------
*/
#if defined(PROBABILISTIC)
	gap[0] = 0.0;
	if (k-j < SCHNITT) {
		for (i=1;i<=k;i++) gap[i] = 0.0;
	} else {
		dum = c[0];
		for (i=1;i<=k;i++) {
			help = 0.0; 
			for (ii=0;ii<i;ii++) help += log(c[ii]);
			help /= (DOUBLE)i;
			help -= (log(pi) + log(2.0)*2.0*(DOUBLE)p/(DOUBLE)i); 
			help += 2.0*loggamma((DOUBLE)i/2.0+1.0)/(DOUBLE)i;
			gap[i] = exp(help);
#if VERBOSE > 0
			printf("%13.10f\n",(double)gap[i]);
#endif
		}
	}
#endif
/*%---------------------------------------------------
*/

/* 
\section{The search loop} 

This may take a very long time. All solutions are enumerated.
*/     
#if defined(PRINTSOLUTION)
	open_solution_file();
#endif
	do {
		loops++;
		if ((stoploops>0)&&(stoploops<=loops)) goto afterloop;
#if VERBOSE > 0              
		if (loops%1000000==0) {
			printf("%ld loops\n",loops);
			fflush(stdout);
		}
#endif
/*   	if ((loops%100000000==0)||(nosolutions>0)) break;*/

/* Calculation of $\tilde{c}_t$. */
		dum = us[t] + y[t];  					
		cs[t] = cs[t+1] + dum*dum*c[t];

/* First test. */

#if defined(PROBABILISTIC)
		if ( (cs[t]<Fd-gap[t]) && (!prune0(fabs(dum),N[t])) )  {   
#else
		if ( (cs[t]<Fd) && (!prune0(fabs(dum),N[t])) )  {   
#endif
			/* Calculation of the new $w_t$. */  
			compute_w(w,bd,dum,t,z);
			if (t>j) {
		
				/* Second test: Test if we are in the affine subspace. */         
				if (pruneprimal(w[t],t,z,Fq,firstlist,firstp)) goto three; 

				/* Third test: Hoelder's inequality. */         
				if ( prune(w[t],cs[t],z,Fq) ) { 
					if (eta[t]==1) goto three;
					eta[t] = 1;
					delta[t] *= -1;
					if (delta[t]*d[t]>=0) delta[t] += d[t];
					us[t] = v[t] + delta[t];
				} else {
					/* We proceed to $t := t-1$. */
					t--;
					eta[t] = 0;
					delta[t] = 0;
					dum = compute_y(mu,us,t,tmax);
					y[t] = dum;
					v[t] = (long)round(-dum);
					us[t] = v[t];
					d[t] = (v[t]>-y[t]) ? -1 : 1; 
				} 
			} else { 
				/* if $t=j$ we can test whether we found a solution: (Third test) */
				if (exacttest(w[j],z,Fq)==1) {          
					/* Output of the solution. */
					onlinesolutiontest(w[t],z,solutionnorm,nom,denom,selection,permutation); 
					if ((stopafter>0)&&(stopafter==nosolutions)) goto afterloop;
				} 
				goto four;
			}
		} else {
three:
			/*  We go back to $t:=t+1$. */
			t++;
			if (tmax<t) tmax = t;
four:		
			if (eta[t]==0) {
				delta[t] *= -1;
				if (delta[t]*d[t]>=0) delta[t] += d[t];
			} else {
				delta[t] += (delta[t]*d[t]>=0) ? d[t] : -d[t] ;
			}
			us[t] = v[t] + delta[t];
		}
	} while (t<=k);	 
/* End of the while-loop. */

afterloop:

#if defined(PRINTSOLUTION)
	close_solution_file();
#endif

/* Here comes some additional output. */
#if VERBOSE > 0              
	printf("Loops: %ld\n",loops);
#endif
	if (((stopafter==nosolutions)&&(stopafter>0))||
		((stoploops<=loops)&&(stoploops>0))) {
		printf("Stopped after number of solutions: %ld\n",nosolutions);
	} else {
		printf("Total number of solutions: %ld\n",nosolutions);
	}	

/* Free the allocated memory. */
	free (us);
	free (cs); 
	free (y);
	free (delta);
	free (d);
	free (first); 
	free (firstlist); 
	free (firstp); 
	free (eta);
	free (v);
	for (l=0;l<s;l++) free (w[l]);
	free (w);
/*%----------------------------------------------------------------
*/
#if defined(PROBABILISTIC)
	free (gap);
#endif
/*%----------------------------------------------------------------
*/
	lllfree(mu,c,N,/*bd,*/s);  
	return 1;
}

/* \section{The pruning functions}
\|exacttest| finally tests if a solution was found. (Maximum norm test)

Normally it should never fail.
*/
int exacttest(DOUBLE *v, int z, DOUBLE Fq) {
	int i;

	if (inequalities!=1) {
		for (i=0;i<z;i++) {
			if (fabs(v[i]) > Fq*(1.0 + EPSILON)) return 0;
		}
		return 1;
	} else {
		for (i=0;i<z-2;i++) {
			if (fabs(v[i]) > Fq*(1.0 + EPSILON)) return 0;
		}
		if (fabs(v[z-2]) > (1.0 + EPSILON)) return 0;
		if (fabs(v[z-1]) > Fq*(1.0 + EPSILON)) return 0;
		return 1;
	}
}

int prune0(DOUBLE li, DOUBLE re) {
	if (li > re) {
		return 1;
	} else {
		return 0;
	}
}

/*  
\|prune| tries to cut the enumeration tree with the help of
 \"Hoelder"'s formula.
*/ 
int prune(DOUBLE *w, DOUBLE cs, int z, DOUBLE Fq) {
       DOUBLE reseite;
       int i;
      
	/* 
	It is tested if
	$$ \tilde{c}_t \leq F_q ||w_t||. $$
	*/
	reseite = -cs*(1.0 - EPSILON)/Fq; 

	for (i=z-1; i>=0; i--) {
		reseite += fabs(w[i]);  
	}
	if (0.0 < reseite) return 0;
	return 1; 
	
}

int pruneprimal(DOUBLE *w, int t, int z, DOUBLE Fq, int *firstlist, int *firstp) {
	int i;
	int f;

	if (inequalities!=1) {
		for (i=0;i<firstlist[firstp[t+1]];i++) {
			f = firstlist[firstp[t+1]+1+i];
			if (fabs(fabs(w[f])-Fq) > EPSILON) return 1;
		} 
		return 0;
	} else {
		for (i=0;i<firstlist[firstp[t+1]];i++) {
			f = firstlist[firstp[t+1]+1+i];
			if ((f<z-2-slackvar)||(f==z-1)) {
				if (fabs(fabs(w[f])-Fq) > EPSILON) return 1;
			} else {
				if (f<z-2) {
					if (fabs(w[f]) > Fq + EPSILON) return 1;
				} else {
					if (fabs(fabs(w[f])-1.0) > EPSILON) return 1;
				}
			}
		} 
		return 0;
	}
}

/* 
\|onlinesolutiontest| is a more explicit test for 
solutions.

It also tests if the resulting design is nontrivial.
*/
int onlinesolutiontest(DOUBLE *w, int z, long c1, long nom, long denom, int selection, int *permutation) {
	int F;
#if defined(PRINTSOLUTION)
	int i,i0,l;
	long u;
	int upper;
#endif

	/* Test if we found a solution. But we are only interested in nontrivial 
   	designs ($\lambda\neq 0$).
	*/
	F = 0;
		
	/* \|labs(b[j][z-1])==1.0| means the solution is inhomogenous.
   	\|b[j][z-2]| indicates a nontrivial solution and \|plus==0| indicates
	   a solution.
	*/   
	if (( fabs(fabs(w[z-1])-(DOUBLE)c1)>EPSILON) || (fabs(w[z-2])<EPSILON)) return 0; 
	/* The next test assures that we only get solutions for the $\lambda$'s
   	which we want to test.
	*/ 
	if (inequalities!=1) {
		if (fabs(fabs(w[z-2])-(DOUBLE)c1)>EPSILON) return 0; 
	} else {
		if (fabs(fabs(w[z-2])-1.0)>EPSILON) return 0; 
	}
	F = 1;
                  
	/* We found a solution, Yipeeeh! */
	if (F==1) {
		nosolutions++;

/* Output of the solution. */           
#if defined(PRINTSOLUTION)
		if (inequalities==1) {
			upper = z-2-slackvar;
		} else {
			upper = z-2;
		}
		/* No preselection of columns */
		if (!selection) {
			for (i=0;i<upper;i++) {
				u=labs((long)round(w[i])-(long)round(w[z-1])*nom)/denom/c1;
				printf("%ld",u);
				fprintf(fp, "%ld",u);
			}  
			if (inequalities==1) {
				for (i=upper;i<z-2;i++) {
					u=labs((long)round(w[i])-(long)round(w[z-1]))/denom;
					printf(" %ld",u);
				}  
			}
		} else {
			i0 = 0;
			for (i=0;i<s_sel-3;i++) {
				u = 0;
				for (l=i0;l<=i;l++) {
					if (permutation[l]==i) {
						u=labs((long)round(w[i0])-(long)round(w[z-1])*nom)/denom/c1;
						i0++;
						break;
					}
				}
				printf("%ld",u);
				fprintf(fp, "%ld",u);
			}  
		}
		printf(" L = %ld\n",(long)labs((long)round(w[z-2])));
		fflush(stdout);
		fprintf(fp, "\n");
		fflush(fp);
#endif
		if (nosolutions%10000==0) {
			printf("%ld\n",nosolutions);
			fflush(stdout);
		}
	}
	return F;
}

/* \section{Orthogonalization}
\subsection{Gram Schmidt orthogonalization}

\|bd| contains the Gram-Schmidt vectors of \|b|.
*/
void gramschmidt(COEFF **b, int k, int z, DOUBLE **mu, DOUBLE **bd, DOUBLE *c, DOUBLE *N, long solutionnorm) {
	int i,l,ii;
	DOUBLE dum;
	
	for (i=0; i<=k; i++) {
		for (l=0;l<z;l++) bd[i][l] = (DOUBLE)b[i][l+1].c;
		N[i] = 0.0;
		for (ii=0;ii<i;ii++) {
			dum = 0.0;
			for(l=0;l<z;l++) dum += (DOUBLE)b[i][l+1].c * bd[ii][l];
			mu[i][ii] = dum / c[ii];
			for (l=0;l<z;l++) bd[i][l] -= mu[i][ii]*bd[ii][l];
		}
	
		c[i] = scalarproductfp(bd[i],bd[i],z);
		for(l=0;l<z;l++) N[i] += fabs(bd[i][l]);
		N[i] /= c[i];
		N[i] *= (DOUBLE)solutionnorm;
#if VERBOSE > 0 	
		printf("%f ",(double)c[i]);
#endif
		N[i] *= 1.0 + EPSILON;
	}
#if VERBOSE > 0 	
	printf("\n\n"); 
#endif
	return;
}

/*
\subsection{Givens rotation}

This is the complete version, i.e. it computes the orthogonal basis,
the Gram Schmidt coefficients, the square of the norm of the orthogonal
basis and the array \|N|.
*/
void givens(COEFF **b, int k, int z, DOUBLE **mu, DOUBLE **bd, DOUBLE *c, DOUBLE *N, long solutionnorm) {
	int i,l,ii;
	int mm;
	DOUBLE d1,d2;
	DOUBLE gc,gs;
	DOUBLE t;

	/* The matrix \|b| is copied to \|mu|.
   	\|bd| is set to a $z\times z$ unity matrix.
	*/   
	for (i=0; i<=k; i++) {
		for (l=0;l<z;l++) {
			mu[i][l] = (DOUBLE)b[i][l+1].c;
		}
	}
	for (i=0; i<z; i++) {
		for (l=0;l<z;l++) {
			bd[i][l] = 0.0;
		}
		bd[i][i] = 1.0;
	}

	/* The Givens rotation */
	for (ii=1; ii<z; ii++) {
		mm = (ii-1 < k) ? ii-1 : k;
		for (i=0;i<=mm;i++) {
			if (mu[i][ii]==0.0) {
				/* Nothing has to be done */
				gc = 1.0;
				gs = 0.0;
			} else {
				/* Stable computation of the
				   rotation coefficients. 
				*/
				if (fabs(mu[i][ii])>=fabs(mu[i][i])) {
					t = mu[i][i]/mu[i][ii];
					gs = 1.0 / SQRT(1.0 + t*t);
					gc = gs * t;
				} else {
					t = mu[i][ii]/mu[i][i];
					gc = 1.0 / SQRT(1.0 + t*t);
					gs = gc * t;
				}
				/* Rotation of \|mu| */
				for (l=i;l<=k;l++) {		
					d1 = mu[l][i];
					d2 = mu[l][ii];
					mu[l][i]  =  gc*d1 + gs*d2;
					mu[l][ii] = -gs*d1 + gc*d2;
				}
				/* Rotation of the matrix $Q^t$ */
				for (l=0;l<z;l++) {		
					d1 = bd[i][l];
					d2 = bd[ii][l];
					bd[i][l]  =  gc*d1 + gs*d2;
					bd[ii][l] = -gs*d1 + gc*d2;
				}
			}
			
		}
	} 

/* Finally some scaling has to be done, since $Q$ is a orthonomal matrix */
	for (i=0;i<=k;i++) {
		c[i] = mu[i][i]*mu[i][i];
		N[i] = 0.0;
		for (ii=0; ii<z; ii++) {
			bd[i][ii] *= mu[i][i];
			N[i] += fabs(bd[i][ii]);
		}
		N[i] /= c[i];
		N[i] *= (DOUBLE)solutionnorm;
		N[i] *= 1.0 + EPSILON;

		for (ii=k; ii>=i; ii--) {
			mu[ii][i] /= mu[i][i];
		}
#if VERBOSE > 0 	
		printf("%f ",(double)c[i]); 
#endif
	}
#if VERBOSE > 0 	
	printf("\n\n"); 
#endif
	return;
}

void givensshort(COEFF **b, int k, int z, DOUBLE **mu, DOUBLE *c) {
	int i,l,ii;
	int mm;
	DOUBLE d1,d2;
	DOUBLE gc,gs;
	DOUBLE t;

	/* The matrix \|b| is copied to \|mu|.
   	\|bd| is set to a $z\times z$ unity matrix.
	*/   
	for (i=0; i<=k; i++) {
		for (l=0;l<z;l++) {
			mu[i][l] = (DOUBLE)b[i][l+1].c;
		}
	}

	/* The Givens rotation */
	for (ii=1; ii<z; ii++) {
		mm = (ii-1 < k) ? ii-1 : k;
		for (i=0;i<=mm;i++) {

			if (mu[i][ii]==0.0) {
				/* Nothing has to be done */
				gc = 1.0;
				gs = 0.0;
			} else {
				/* Stable computation of the
				   rotation coefficients. 
				*/
				if (fabs(mu[i][ii])>=fabs(mu[i][i])) {
					t = mu[i][i]/mu[i][ii];
					gs = 1.0 / SQRT(1.0 + t*t);
					gc = gs * t;
				} else {
					t = mu[i][ii]/mu[i][i];
					gc = 1.0 / SQRT(1.0 + t*t);
					gs = gc * t;
				}
				/* Rotation of \|mu| */
				for (l=i;l<=k;l++) {		
					d1 = mu[l][i];
					d2 = mu[l][ii];
					mu[l][i]  =  gc*d1 + gs*d2;
					mu[l][ii] = -gs*d1 + gc*d2;
				}
			}
			
		}
	} 

/* Finally some scaling has to be done, since $Q$ is a orthonomal matrix */
	for (i=0;i<=k;i++) {
		c[i] = mu[i][i]*mu[i][i];
		for (ii=k; ii>=i; ii--) {
			mu[ii][i] /= mu[i][i];
		}
	}
	return;
}

/* \section{The $\Gamma$ Function}
The logarithm of the $\Gamma$-Function computed with the
Euler-McLaurin summation formular. */

DOUBLE laurin(DOUBLE x) {
#define K1  0.918938533204672741780329
#define K2  0.083333333333333333333333
#define K3  0.002777777777777777777777
#define K4  0.000793650793650793650793
#define K5  0.000595238095238095238095
#define K6  0.000841750841750841750841
#define K7  0.001917526917526917526917
#define K8  0.006410256410256410256410
#define K9  0.029550652951002120971679
#define K10 0.179681226611137390136718
#define K11 1.39243221282958984375
#define MM  13

	DOUBLE y;
	
	if (x<=0.0) {
		return 0.0;
	}
	y = 1.0/(x*x);
	x = (x-0.5)*log(x) - x + K1 + (K2-y*(K3-y*(K4-y*(K5-y*(K6-y*(K7-y*(K8-y*(K9-y*(K10-y*K11)))))))))/x;
	return x;
}

/* The logarithm of the $\Gamma$ function is computed. */
DOUBLE loggamma(DOUBLE x) {
	int n;
	DOUBLE y;
	int i;
	
	if ((x<=0.0) || (x>1.0e70)) {
		return 0.0;
	}
	if (x>1.0e10) {
		return x*(log(x)-1.0);
	}
	if (x>=MM) {
		return laurin(x);
	}
	n = MM - (long)(x);
	
	/* After some tests \|laurin| is called. */
	y = x - floor(x) + MM;
	y = laurin(y);
	for (i=0;i<n;i++) {
		y -= log(x+i);
	}
	return y;
}

/* end of file \|enumall.c| */

