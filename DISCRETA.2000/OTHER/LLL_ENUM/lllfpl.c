/*
% $Header: /usr/local/cvsroot/designs/lllfpl.c,v 1.6 1997/11/11 19:17:41 alfred Exp $
 \input cnoweb
 \title{LLL-reduction using floating point arithmetic}
 \job{lllfpl.c}
 \synopsis{LLL-algorithm using floating point arithmetic 
 with correction of rounding errors, see {\it Schnorr}, {\it Euchner}, 
 'Fundamentals of Computation Theory' (1991).

History:
\item{--} Written by Alfred Wassermann, 24.12.94.
\item{--} Pruned Gauss Enumeration and Sparse structure added 6.6.1995.
\item{--} Additional parameter for pruned Gauss enumeration ($p$).
\item{--} 23.1.1996: Bug found and corrected in Step 3.
\item{--} 28.10.1996: Minor corrections and removing of old and unused code.
\item{--} 14.3.1997: scalarproductl(fp) now returns double-value.

 Input: 
 \item{} \|b|,\|s|,\|z|: \|COEFF|-matrix with \|s| columns and \|z| rows.
 The columns of \|b| build the basis of the lattice.
 \item{} \|mu|, \|c|: Return the Gram Schmidt coefficients and the 
   norms of the orthogonal basis vectors.
 \item{} \|delta|: determines the reduction quality.
 Originally $\delta = 3/4$. Schnorr takes $\delta= 0.99$.

 Compiler options:
 \item{} \|-DDEEPINSERT|: Step 4 is replaced by ''deep insertions'', 
	 see {\it Schnorr}.

In \|llll.h| several other header files are loaded.
The algorithm uses \|long|-integers. If there are errors, the
multiprecision version has to be used.
The Unix-libraries \|math| and \|malloc| are used.}
\section{Main program: subsetsum}
\section{Introduction}  
The standard libraries, the header file with the return-values are included.

If \|DOUBLE| is defined as type \|double| then
the \|DOUBLE| constant \|TWOTAUHALF| is set to $2^{52/2}$, where the machine
precision of a double variable is 53 bit for the 486.
If \|DOUBLE| is defined as type \|long double| then
the \|DOUBLE| constant \|TWOTAUHALF| is set to $2^{64/2}$, where the machine
precision of a long double variable is 64 bit for the 486.
If \|DOUBLE| is defined as type \|long double| then
the \|DOUBLE| constant \|TWOTAUHALF| is set to $2^{128/2}$, where the machine
precision of a long double variable is 128 bit for the CRAY.

{\bf Attention:} the first index of \|b[i][j]| are the columns of
$b_{j,i}$, the second index is the row number!!! 
*/
#include "llll.h"
#define TWOTAUHALF 67108864.0 
/* *< #define TWOTAUHALF 4294967296.0  >* */
/* *< #define TWOTAUHALF 18446744073709551616.0 >* */

/*
Now we begin with the LLL algorithm with floating point corrections.
*/
int lllfp (COEFF **b, DOUBLE **mu, DOUBLE *c, DOUBLE *N, /* DOUBLE **bs, */ int start, int s, int z, DOUBLE delta, int KMz, long v1, long v2, long c1, long nom, long denom, int flag)
{
/* Additional input variables:
   \item{} \|N|: the $l_2$-norm of $b_i$ (rounded).
   \item{} \|c|: the square of the $l_2$-norm of $\hat{b}_i$, where
		 $\hat{b}_i$ are the Gram-Schmidt-vectors.
   \item{} \|mu|: $\mu$, the lower triangle matrix, containing the
		 Gram-Schmidt-coefficients of the orthogonalization.
   \item{} \|bs|: a matrix containing the rounded vectors $b_i$.
                  This matrix is not used in this version. It is 
                  only useful in the multiprecision version.
*/

	int      i,j,k;
	int      Fc, Fr, Fi;
	DOUBLE   mus;
	long     musvl;
	COEFF    *swapvl;
/* \|swapd| is only needed if \|bs| is needed. */
	/* *< DOUBLE   *swapd; >* */
	DOUBLE   ss;
	DOUBLE   cc;
	int      counter;

	DOUBLE	*muc;
	int 	ii,iii;
	COEFF	*bb;
#if defined(DEEPERINSERT)
	long	uu[3];
	DOUBLE	NN;
#endif

/*
Are the dimensions correct?
*/

	if ((z <= 1) || (s <= 1)) return(Matrix_zu_klein);
	
/* New!!! \|muc| is used to save some multiplications in Step 2.*/
	muc=(DOUBLE*)calloc(s,sizeof(DOUBLE));
	if (muc == 0) return(zu_wenig_Speicherplatz);

/*
\section{Step 1}
We begin with stage $k=1$. Schnorr's original algorithm runs
from $k=2$ to $k=s$. Here we have $k=1,\ldots,s-1$.

Initially the $b_i$'s are rounded and then their norms,
$N_i = ||b'_i||$, are computed.

Then we do steps 2 -- 4 while  $k\leq s$.
*/
	k = start;
	if (k<1) k = 1;
	for (i=k-1;i<s;++i) {
		ss = 0.0;
/* In the \|long| version \|bs| is'nt really needed. */
/* *< 		for (j=0;j<z;++j) {
			bs[i][j] = (DOUBLE)b[i][j+1].c;
			ss += bs[i][j]*bs[i][j];
                }
>* */
/* *< 		ss = (DOUBLE)norml(b[i],z);               >* */ 
		ss = normfp(b[i]/*,z*/);                
		N[i] = SQRT(ss);
	}

	/*
	   \|counter| is just kept for statistical purposes.
	*/
	counter = 0;
	while (k<s) {
#if VERBOSE > 2
			if ((counter%500)==0) {printf("LLL: %d \n",counter);} 
			counter++;
#endif

/*
\section{Step 2}
computation of $\mu_{k,1},\ldots,\mu_{k,k-1}$ and $c_k = ||\hat{b}_k'||^2$.
*/

		if (k==1) c[0] = N[0]*N[0];
		c[k] = N[k]*N[k];
		for (j=0;j<k;j++) {
/* It seems that on my 486 PC \|scalarproductfp| takes twice the time of
   \|scalarproductl|. */
/* *<   			ss = scalarproductfp(bs[k],bs[j],z);
			if (fabs(ss)< N[k]*N[j]/TWOTAUHALF) {
				ss = (DOUBLE)scalarproductl(b[k],b[j],z);
			} 
>* */
			
       			ss = scalarproductlfp(b[k],b[j]/*,z*/); 
       			
/* The original commands. */       			
/* *<			for (i=0;i<j;i++) ss -= mu[j][i]*mu[k][i]*c[i]; >* */
/* Here some multiplications are saved. */			
			for (i=0;i<j;i++) ss -= mu[j][i]*muc[i]; 
			
/* The next 5 lines are just for test purposes. They are in comments because
   so far I had no problems. */
/* *<			if (fabs(c[j])<1.0e-10) {
				printf("c[j] zu klein-Ende\n");
				return 0;
			}
>* */			
/* \|muc| is updated. */			
			muc[j] = ss; 
			mu[k][j] = ss/c[j];
			c[k] -= ss*mu[k][j];
		}  /* for $j$ */

/*
\section{Step 3}
size-reduction of $b_k$.
*/
		Fc = 0; Fr = 0;                 /* $F_c = F_r = {\rm false}$  */
		for(j=k-1;j>=0;j--) {
			if (fabs(mu[k][j])>0.5) {
				mus = round(mu[k][j]);
				musvl = (long)mus;
				Fr = 1;           /* $F_r = {\rm true}$: 
				                     Reduction has to be done. 
				                  */
				if (fabs(mus)>TWOTAUHALF) {
					 Fc = 1;  /* $F_c = {\rm true}$ 
					             We will step back due to
					             rounding errors.
					          */
				}
/* Update $b_k$ with $b_k = b_k - \mu b_j$ and
   $\mu_{k,i}$ with $\mu_{k,i} = \mu_{j,i}\lceil\mu_{k,j}\rfloor$.
   This is done using the sparse structure of the vectors.
   
   In order to save multiplications, we treat the cases where
   $\lceil\mu_{k,j}\rfloor = \pm 1$ separately.
*/
				switch (musvl) {
				/* $\lceil\mu_{k,j}\rfloor = 1$: */
				case 1:
					/* The original loop: */
					/* *< for(i=1;i<=z;i++) b[k][i].c -=  b[j][i].c; >* */
					i=b[j][0].p;
					if (i!=0) {
						do {
							bb=&(b[k][i]);
							bb->c -=  b[j][i].c;
							iii = bb->p;
							if ((b[k][i-1].p!=i)&&(bb->c!=0)) {
								for (ii=i-1;(b[k][ii].p==iii)&&(ii>=0);ii--) b[k][ii].p = i;
							}
							else if (bb->c==0) {
								for (ii=i-1;(b[k][ii].p==i)&&(ii>=0);ii--) b[k][ii].p = iii;
							}
							i = b[j][i].p;
						} while (i!=0);
					}
					for(i=0;i<j;i++) mu[k][i] -= mu[j][i];
					break;
					
				/* $\lceil\mu_{k,j}\rfloor = -1$: */
				case -1:
					/* The original loop: */
					/* *< for(i=1;i<=z;i++) b[k][i].c +=  b[j][i].c; >* */
					i=b[j][0].p;
					if (i!=0) {
						do {
							bb=&(b[k][i]);
							bb->c +=  b[j][i].c;
							iii = bb->p;
							if ((b[k][i-1].p!=i)&&(bb->c!=0)) {
								for (ii=i-1;(b[k][ii].p==iii)&&(ii>=0);ii--) b[k][ii].p = i;
							}
							else if (bb->c==0) {
								for (ii=i-1;(b[k][ii].p==i)&&(ii>=0);ii--) b[k][ii].p = iii;
							}
							i = b[j][i].p;
						} while (i!=0);
					}
					for(i=0;i<j;i++) mu[k][i] += mu[j][i];
					break;
					
				/* $\lceil\mu_{k,j}\rfloor \neq \pm 1$: */
				default:
					/* The original loop: */
					/* *< for(i=1;i<=z;i++) b[k][i].c -=  b[j][i].c*musvl; >* */
					i=b[j][0].p;
					if (i!=0) {
						do {
							bb=&(b[k][i]);
							bb->c -=  b[j][i].c*musvl;
							iii = bb->p;
							if ((b[k][i-1].p!=i)&&(bb->c!=0)) {
								for (ii=i-1;(b[k][ii].p==iii)&&(ii>=0);ii--) b[k][ii].p = i;
							}
							else if (bb->c==0) {
								for (ii=i-1;(b[k][ii].p==i)&&(ii>=0);ii--) b[k][ii].p = iii;
							}
							i = b[j][i].p;
						} while (i!=0);
					}
					for(i=0;i<j;i++) mu[k][i] -= mu[j][i]*mus;
				}
				
				mu[k][j] -= mus;
				
				/* Test for a solution. */
  				solutiontest(b[k], z, KMz, v1, v2, c1, nom, denom,flag);
			}   /* end if $|\mu_{k,j}|$ */
		}  /* end for $j$ */

		if (Fr==1) {
/* With \|long|-arithmetics we don`t use \|bs|. */
/* *<		
			N[k] = 0.0;
			for(i=0;i<z;i++) {
				bs[k][i] = (DOUBLE)b[k][i+1].c;
				N[k] += bs[k][i]*bs[k][i];
			} >* */   /* for $i$ */

		/* *< 	N[k] = (DOUBLE)norml(b[k],z);                >* */
			N[k] = normfp(b[k]/*,z*/);                
			N[k] = SQRT(fabs(N[k]));
		}       /* if $F_r$ */

/* Now the correction due to rounding errors. If $F_c$ is
true, then the stage is decreased by one. Otherwise
the new value of $b'_k$ is computed and we proceed to step 4.
*/

		if (Fc==1) k = (k-1 > 1) ? k-1 : 1;    /* $k=\max(k-1,1)$ */
		else {

/* At the end of step 3 we test if $b_k$ is linear dependent.
   Then if $||N_k||<{1\over 2}$ we shift $b_k$ to the last column of the
   matrix and restart \|lllfp| with $s:= s-1$.
*/
			if (fabs(N[k])<0.5) {
				   swapvl = b[k];
				   ss = N[k];
				   /* *< swapd = bs[k]; >* */
				   for(i=k+1;i<s;i++) {
				   	b[i-1] = b[i];
				   	N[i-1] = N[i];
				   /* *<	bs[i-1] = bs[i]; >* */
				   }
				   b[s-1] = swapvl;
				   N[s-1] = ss;
				   /* *< bs[s-1] = swapd; >* */
		
				   /* Now we decrease the size of the
				      lattice basis. */
				   s = s-1;
				   /* With setting $k := 1$ and \|continue| we
				      restart \|lllfp|.
				   */
				   k = 1;
				   continue;
				}

/*
\section{Step 4}
\subsection{New Step 4}
Will be used if the compiler flag \|-DDEEPINSERT| is set.
The vector $b_k$ will be inserted left as possible.

Swapping of the matrix columns is done entirely with swapping of
pointers. Hope it works. \|;-)|
*/
#if defined(DEEPINSERT)

			cc = N[k]*N[k];
			j = 0;
			Fi = 0;
			while (j<k) {
				if (delta*c[j]<=cc) {
					cc -= mu[k][j]*mu[k][j]*c[j];
					j++;
				}
				else {
/* insert $b_k$ at position $j$,
   insert $b'_k$ at position $j$,
   insert $N_k$ at position $j$: */
					swapvl = b[k];
					/* *< swapd = bs[k]; >* */
					ss = N[k];
					for(i=k-1;i>=j;i--) {
						b[i+1] = b[i];
						/* *< bs[i+1] = bs[i]; >* */
						N[i+1] = N[i];
					}						
					b[j] = swapvl;
					/* *< bs[j] = swapd; >* */
					N[j] = ss;

					Fi = 1;       /* $F_i= {\rm true}$ */
					break;
				}     /* end if $\delta\ldots$ */
			}             /* end while  */
			if (Fi==1) k = (j > 1) ? j : 1;
			else k++;
#else

/*
\subsection{Original Step 4}
Either swap the vectors $b_{k-1}$ and  $b_k$ or increment $k$.
Swapping of the matrix columns is done entirely with swapping of
pointers. 
*/

			if (delta*c[k-1] > c[k]+mu[k][k-1]*mu[k][k-1]*c[k-1]) {
/* swap $b_k, b_{k-1}$: */
				swapvl = b[k];
				b[k] = b[k-1];
				b[k-1] = swapvl;
/* swap $b'_k, b'_{k-1}$: */
				swapd = bs[k];
				bs[k] = bs[k-1];
				bs[k-1] = swapd;
/* swap $N_k, N_{k-1}$: */
				ss = N[k];
				N[k] = N[k-1];
				N[k-1] = ss;
/* $k = \max(k-1,1)$ */
				k = (k-1 > 1) ? k-1 : 1;
			}       /* end if $\delta\ldots$ */
			else k++;
#endif  /* \|DEEPINSERT|  */
		}   /* else $F_c$ */
/* end of loop */
	}
	/* while $k<s$ */

#if VERBOSE > 1
		if (counter>1) {
			printf("Schleifen: %d\n", counter);
		}
#endif	

	free(muc);
	return(regulaerer_Durchlauf);
}

/*
\section{Additional subroutines}
\subsection{Memory allocation}
*/
int lllalloc(DOUBLE ***mu, DOUBLE **c, DOUBLE **N, /*DOUBLE ***bs,*/ int s, int z)
{
	int i,j;

/* {\bf Allocation}
 For all matrices except $\mu$ we have to obey the following rule:
The first index gives the column number, the second index the row number.

Are the dimensions correct?
*/

	if ((z < 1) || (s < 1)) return(Matrix_zu_klein);

/* Allocation for c: \|s| \|DOUBLE|-entries */
	(*c)=(DOUBLE*)calloc(s,sizeof(DOUBLE));
	if ((*c) == 0) return(zu_wenig_Speicherplatz);
/* Allocation for N: \|s| \|DOUBLE|-entries */
	(*N)=(DOUBLE*)calloc(s,sizeof(DOUBLE));
	if ((*N) == 0) return(zu_wenig_Speicherplatz);
/* Allocation for mu.  2. Variants:
\item{--} triangular matrix for the original Gram Schmidt orthogonalization,
\item{--} full matrix for Givens rotation. Of course also Givens rotation
          can be done with a triangular matrix. But since memory is not the
          problem I'm somewhat lazy.

First the full matrix for Givens rotation. 
*/

#if defined(GIVENS)
	(*mu)=(DOUBLE**)calloc(s,sizeof(DOUBLE*));
	if ((*mu) == 0)
		{
		free((*c));
		free((*N));
		return(zu_wenig_Speicherplatz);
		}
	for(i=0;i<s;++i)
		{
		(*mu)[i]=(DOUBLE*)calloc(z,sizeof(DOUBLE));
		if ((*mu)[i] == 0)
			{
			for (j=1;j<i;++j) free((*mu)[j]);
			free((*mu));
			free((*N));
			free((*c));
			return(zu_wenig_Speicherplatz);
			}
		}
#else		
/* Now the Gram Schmidt version.*/
	(*mu)=(DOUBLE**)calloc(s,sizeof(DOUBLE*));
	if ((*mu) == 0)
		{
		free((*c));
		free((*N));
		return(zu_wenig_Speicherplatz);
		}
	for(i=1;i<s;++i)
		{
		(*mu)[i]=(DOUBLE*)calloc(i,sizeof(DOUBLE));
		if ((*mu)[i] == 0)
			{
			for (j=1;j<i;++j) free((*mu)[j]);
			free((*mu));
			free((*N));
			free((*c));
			return(zu_wenig_Speicherplatz);
			}
		}
#endif
		
/* Allocation for \|bs| */
/* \|bs| is not needed in the \|long|-version. */
/* *<
	(*bs)=(DOUBLE**)calloc(s,sizeof(DOUBLE*));
	if ((*bs) == 0)
		{
		for (i=1;i<s;++i) free((*mu)[i]);
		free((*mu));
		free((*c));
		free((*N));
		return(zu_wenig_Speicherplatz);
		}
>* */		
/* *<
	for(i=0;i<s;++i)
		{
		(*bs)[i]=(DOUBLE*)calloc(z,sizeof(DOUBLE));
		if ((*bs)[i] == 0)
			{
			for (j=0;j<i;++j) free((*bs)[j]);
			free((*bs));
			for (j=1;j<s;++j) free((*mu)[j]);
			free((*mu));
			free((*c));
			free((*N));
			return(zu_wenig_Speicherplatz);
			}
		}
>* */		
	return(regulaerer_Durchlauf);
}
/*
\subsection{Free the allocated memory}
*/
int lllfree(DOUBLE **mu, DOUBLE *c, DOUBLE *N, /*DOUBLE **bs,*/ int s)
{
	int i;

	for (i=1;i<s;++i) free(mu[i]);
	free(mu);
	free(c);

/* \|bs| is not needed in the \|long|-version. */
/* *<
	for (i=0;i<s;++i) free(bs[i]);
	free(bs);
>* */	
	free(N);
	return(regulaerer_Durchlauf);
}

/*
\subsection{Scalarproducts}
Two scalarproducts are defined: One for floatingpoint arithmetic, the
other for \|long|'s.
*/

DOUBLE scalarproductfp (DOUBLE *v, DOUBLE *w , int n)
{
	int i;
	DOUBLE r;

	r = 0.0;
	for (i=n-1;i>=0;i--) r += v[i]*w[i];
	return (r);
}

/* The integer scalarproduct is able to use the sparse structure of the
   vectors. But it strongly depends on the hardware and the compiler.
   
   It seems that on the 486 PC the \|long| variant is twice as fast as
the \|DOUBLE| variant.

This is the new variant which uses the sparse structure of the input 
vectors.
*/
long scalarproductl (COEFF *v, COEFF *w /*, int n*/)
{
	long	erg;
	long	t1,t2;
	COEFF	*vv,*ww;

	erg = 0;
	/* The old instruction. */
	/* *< for (i=n-1;i>=0;i--)  erg += v[i+1].c*w[i+1].c; >* */
	
	t1 = v[0].p;
	t2 = w[0].p;
	if ((t1==0)||(t2==0)) return 0;

/* It depends on the compiler and the machine whether this loop
is fast enough to gain speed against the plain scalarproduct evaluation.
$t_1$ runs through the vector $v$, $t_2$ runs through the vector $w$.
*/	
	do {
		/* If $t_1 \neq t_2$ we don't have to work. 
		   We can step forward. 
		*/
		if (t2>t1) {
			t1 = v[t2-1].p;
			if (t2!=t1) { 
				if (t1==0) break;
				t2 = w[t2].p; 
				if (t2==0) break;
			}
			else goto gleich;
		}
		else if (t2<t1) {
			t2 = w[t1-1].p;
			if (t2!=t1) { 
				if (t2==0) break;
				t1 = v[t1].p; 
				if (t1==0) break;
			}
			else goto gleich;
		}
		else {
/* Only in the case $t_1=t_2$ we have to multiply. */
gleich:
			vv = &(v[t1]);
			ww = &(w[t2]);
	     		erg += vv->c * ww->c;
			t1 = vv->p;
			if (t1==0) break;
			t2 = ww->p;
			if (t2==0) break;
		}
	}
	while (1);

	return (erg);
}

DOUBLE scalarproductlfp (COEFF *v, COEFF *w /*, int n*/)
{
	DOUBLE	erg;
	long	t1,t2;
	COEFF	*vv,*ww;

	erg = 0.0;
	/* The old instruction. */
	/* *< for (i=n-1;i>=0;i--)  erg += (DOUBLE)(v[i+1].c*w[i+1].c); >* */
	
	t1 = v[0].p;
	t2 = w[0].p;
	if ((t1==0)||(t2==0)) return 0;

/* It depends on the compiler and the machine whether this loop
is fast enough to gain speed against the plain scalarproduct evaluation.
$t_1$ runs through the vector $v$, $t_2$ runs through the vector $w$.
*/	
	do {
		/* If $t_1 \neq t_2$ we don't have to work. 
		   We can step forward. 
		*/
		if (t2>t1) {
			t1 = v[t2-1].p;
			if (t2!=t1) { 
				if (t1==0) break;
				t2 = w[t2].p; 
				if (t2==0) break;
			}
			else goto gleich;
		}
		else if (t2<t1) {
			t2 = w[t1-1].p;
			if (t2!=t1) { 
				if (t2==0) break;
				t1 = v[t1].p; 
				if (t1==0) break;
			}
			else goto gleich;
		}
		else {
/* Only in the case $t_1=t_2$ we have to multiply. */
gleich:
			vv = &(v[t1]);
			ww = &(w[t2]);
	     		erg += (DOUBLE)(vv->c) * (DOUBLE)(ww->c);
			t1 = vv->p;
			if (t1==0) break;
			t2 = ww->p;
			if (t2==0) break;
		}
	}
	while (1);

	return (erg);
}

/* Computation of the square of the norm. This is also done using
   the sparse structure. We have two (sparse) versions: 
   \|norml| computes an
   \|integer| result and \|normfp|  a \|DOUBLE| result. Both procedures
   have an \|integer| vector as input.
*/

long norml (COEFF *v/*, int n*/)
{
	long	erg;
	long	t;
	COEFF	*vv;

	erg = 0;
	
	t = v[0].p;
	if (t==0) return 0;

	do {
		vv = &(v[t]);
     		erg += vv->c * vv->c;
		t = vv->p;
	}
	while (t!=0);

	return (erg);
}

DOUBLE normfp (COEFF *v/*, int n*/)
{
	DOUBLE	erg;
	long	t;
	COEFF	*vv;

	erg = 0.0;
	
	t = v[0].p;
	if (t==0) return 0;

	do {
		vv = &(v[t]);
     		erg += (DOUBLE)(vv->c) * (DOUBLE)(vv->c);
		t = vv->p;
	}
	while (t!=0);

	return (erg);
}

/* \subsection{Init procedure for \|COEFF|-vectors.}

    The last value of $r$: \|v[0].c| is thrown away.
*/
void coeffinit(COEFF *v, int z)
{
	short r = 0;
	short i;
	for (i=z;i>=0;i--) {
		v[i].p = r;
		if (v[i].c!=0) r=i;
	}
	return;
}

/*
\subsection{rounding procedure}
\|round(DOUBLE r)| denotes $\lceil r - 1/2\rceil$ with
\|round(r)|$ = r-1/2$ for half integers.
*/

DOUBLE round(DOUBLE r)
{
 return( ceil(r-0.5) );
}

/*
\subsection{Abstract LLL}
This routine calls \|lllfp| with $\delta=$\|LLLCONST_LOW|.
It is meant for generic basis reduction without nearestpoint
calculation afterwards.
*/
int lll (COEFF **b, int s, int z,int KMz, long v1, long v2, long c1, long nom, long denom, int flag)
{
	DOUBLE **mu;
	DOUBLE *c;
	DOUBLE *N;
/*	DOUBLE **bs;*/
	DOUBLE delta = LLLCONST_LOW;
	int r;

	lllalloc(&mu,&c,&N,/*&bs,*/s,z);
	r = lllfp(b,mu,c,N,/*bs,*/1,s,z,delta,KMz,v1,v2,c1,nom,denom,flag);
	lllfree(mu,c,N,/*bs,*/s);

	return(r);
}

/*
\section{Input/Output routines}

 \|ausgabe| writes the matrix \|b|.
 \|b| is a $z\times s$-matrix.

{\bf Attention:} the first index of $b_{i,j}$ or \|b[i][j]| is the
column number, the second index is the row number!!!
*/

void ausgabe(COEFF **b, int s, int z)
	{
	int i,j;

	printf("  Spaltenzahl: %d\n",s);
	printf("  Zeilenzahl:  %d\n",z);
	for (j=0;j<s;++j) {
        	for (i=1;i<=z;++i) printf("%ld ",b[j][i].c);
		printf("\n");
	}
	printf("\n");
	fflush(stdout);
	return;
	}

/*
 \|vectorausgabe| writes the vector \|v|.
 \|v| is a $z$-vector.
*/

void vectorausgabe(COEFF *v, int z)
	{
	int i;

	printf("  Zeilenzahl:  %d\n",z);
	for (i=1;i<=z;++i) printf("%ld ",v[i].c);
        printf("\n");
        fflush(stdout);
	return;
	}

/* \subsection{Output and input intermediate results}

    \|putlattice| write a complete lattice basis to the file 
    \|lll.output|.
    The first line contains the number of rows and the number
    of colums of the lattice. Moreover it contains the numer of rows and columns
    of the underlying Kramer Mesner matrix.
    The last column is not writte because it is equal to zero.
*/

void putlattice(COEFF **b, int s, int z, int KMz, int KMs)
{
	FILE *out_txt;
	int i,j;
	out_txt = fopen("lll.output", "w");
	fprintf(out_txt,"%d %d %d %d\n",z,s,KMz,KMs);
	for (j=1;j<=z;j++) {
		for (i=0;i<s-1;i++) fprintf(out_txt,"%ld\t",b[i][j].c);
		fprintf(out_txt,"\n");
	}			
	fclose(out_txt);
	return;
}	

/* \|getlattice| reads a complete lattice basis as it is written by \|putlattice|.
   There are two versions:
   \item{--} \|flag|${} = 0$: It is completely read.
   \item{--} \|flag|${} = 1$: The columns which contain
    nonzero entries in the rows belonging to the underlying Kramer Mesner matrix
    are deleted. Then the rows belonging to the Kramer Mesner matrix are zero.
    These are deleted too.
    
    The column number $s-1$ is not read it has to be set to zero.
*/    
void getlattice(COEFF **b, int *s, int *z, int flag) {
	FILE *in_txt;
	int i,j,KMz,KMs;
	long test;
	in_txt = fopen("lll.output", "r");
	
/* This version assumes that \|b| is already allocated. */
	for (j=0;j<=(*z);j++) 
		for (i=0;i<(*s);i++) b[i][j].c = 0;
		
/* First the lattice basis is completely read. */		
	fscanf(in_txt,"%d%d%d%d\n",&(*z),&(*s),&KMz,&KMs);
	
	for (j=1;j<=(*z);j++) {
		for (i=0;i<(*s)-1;i++) fscanf(in_txt,"%ld",&(b[i][j].c));
	}
	fclose(in_txt);
		
/* Here comes the deleting part: */
	if (flag==1) {
		j=0;
		do {
			test = 0;
			for(i=1;i<=KMz;i++) test += labs(b[j][i].c);
			if ((test==0) /* *< &&(labs(b[j][(*z)].c)<=1) >* */ ) j++;
			else {
				for (i=j+1;i<(*s);i++) b[i-1] = b[i];
				(*s)--;
			}
		} while (j<(*s)-1);
		
		/* Test if there are solutions possible. */
		test = 0;
		for (i=0;i<(*s);i++) if (b[i][(*z)].c!=0) {
			test = 1;
			break;
		}
		if (test==0) {
			printf("No solutions after phase 1.\n");
			exit(0);
		}			

/* Now the rows are deleted. */
		for (j=0;j<(*s);j++)
			for (i=KMz+1;i<=(*z);i++) b[j][i-KMz].c = b[j][i].c;
		(*z) -= KMz;
	}
	
	/* Update of the sparse structure. */
 	for (j=0;j<(*s);j++) coeffinit(b[j],(*z));
	
	return;
}

/* Delete unnecessary columns and rows after the first phase
of the algorithm without output on file. */

void cutlattice(COEFF **b, int *s, int *z, int KMz) {
	int j,i;
	long test;
	
	j=0;
	do {
		test = 0;
		for(i=1;i<=KMz;i++) test += labs(b[j][i].c);
		if (test==0) j++;
		else {
			for (i=j+1;i<(*s);i++) b[i-1] = b[i];
			(*s)--;
		}
	} while (j<(*s)-1);
		
	/* Test if there are solutions possible. */
	test = 0;
	for (i=0;i<(*s);i++) if (b[i][(*z)].c!=0) {
		test = 1;
		break;
	}
	if (test==0) {
		printf("No solutions after phase 1.\n");
		exit(0);
	}			

	/* Now the rows are deleted. */
	for (j=0;j<(*s);j++)
		for (i=KMz+1;i<=(*z);i++) b[j][i-KMz].c = b[j][i].c;
	(*z) -= KMz;

	/* Update of the sparse structure. */
 	for (j=0;j<(*s);j++) coeffinit(b[j],(*z));
	
	return;
}

/* \subsection{analyse lattice size}

\|analyse| computes and outputs the norms of the lattice vectors and outputs
the sum of the norms.
Returns the sum of the norms.
*/

double analyse(COEFF **b, int s, int z)
	{
	double ss,r,x;
	int i,j;

	ss = 0.0;
	for(j=0;j<s;j++) {
		r = 0.0;
		for (i=1;i<=z;i++) {
			x = (double)b[j][i].c;
			r += x*x;
		}
		/* *<  printf("%e ",r); >* */
		ss += r;
	}
#if VERBOSE > 1
	printf("Analyse: %13.10f\n",ss);
#endif	
	return(ss);
	}

/*
\section{\"bkz": Blockwise Korkine Zolotarev Reduction}
There are 2 parameters $\beta$ and $\delta$ with ${1\over 4}<\delta<1$ (same
as in LLL) and $2<\beta<m$ (the blocksize).
It uses the subroutines \|lllfp| and \|enumerate|.

$p$ controls the pruning in the pruned Gauss enumeration.
*/

int bkz (COEFF **b, int s, int z, DOUBLE delta, int beta, int KMz, long v1, long v2, long c1, long nom, long denom, int flag, int p)
{

	DOUBLE **mu;
	DOUBLE *c;
	DOUBLE *N;
/*	DOUBLE **bs; */

/*
   $j$ is cyclically shifted through $0,1,\ldots,s-2$.
   \|zaehler| counts the number of positions $j$ that satisfy
   $\delta||\hat{b}_j||^2 \leq \lambda_1(\pi_j(L(b_j,\ldots,b_k)))$.
   \|zaehler| is reset to $0$ if the inequality does not hold for $j$.
*/
     int zaehler;
     int j,h,k;
     int i,l;
     DOUBLE cjdach;
     long *u;
     COEFF *swapvl;
     
     if (s<2) {
        printf("The LLL algorithm didn't succeed in computing a basis of the kernel.\n");
        printf("Please increase c0 (the first parameter)!\n");
	return 0;
     }

/* allocate memory for \|u| */
	u = (long*)calloc(s+1,sizeof(long));
	for (i=0;i<s+1;i++) u[i] = 0;

	lllalloc(&mu,&c,&N,/*&bs,*/s,z);
/*
 \"Step 1".
 $\mu$ and $c$ have to be allocated with $s+1$ columns.
*/
#if VERBOSE > 2
     	printf("Start of bkz.\n");
#endif
     lllfp(b,mu,c,N,/*bs,*/1,s,z,delta,KMz,v1,v2,c1,nom,denom,flag);

     zaehler = 0;
     j = 0;

   k = (j+beta-1 < s) ? j+beta-1 : s;
   if (j<k) {
     while (zaehler<s-1) {
/* Debugging output. */
#if VERBOSE > 2
		printf("bkz: %d \n",zaehler);  
#endif

	j += 1;
	/*  $k := \min(j+\beta-1,s)$  */
	k = (j+beta-1 < s) ? j+beta-1 : s;
	if (j==s) {
		j = 1;
		k = (beta < s) ? beta : s;
	}

/*
  \"Step 2"
*/
       cjdach = enumerate(mu,c,u,s,j-1,k-1,p);

/*
  \"Step 3"
*/
	/*  $h := \min(k+1,s)$  */
	h = (k+1 < s) ? k+1 : s;
	if (fabs(delta*c[j-1]) > fabs(cjdach)) {
/* 
   Now we found a new vector whose corresponding Gram Schmidt vector
   is shorter than $c_j$ by at least a factor $\delta$.
*/
#if VERBOSE > 1
		printf("down %d  %d \n",j,k);  
#endif
		swapvl = b[s];
		for (i=s-1;i>=j-1;i--) b[i+1] = b[i];
		b[j-1] = swapvl;

/* Write the new linear combination at position $j-1$ after shifting all 
   vectors behind $j-1$ one position to the right.
*/   
		for (l=1;l<=z;l++) b[j-1][l].c = 0;
		for (i=j;i<=k;i++) {
			for (l=1;l<=z;l++) b[j-1][l].c += b[i][l].c*u[i-1];
		}
		/* Lazy update of the coefficients of the new $b_{j-1}$. */
		coeffinit(b[j-1],z);
		
#if VERBOSE > 2
			for(i=j-1;i<k;i++) printf("%ld ",u[i]);
			printf("\n"); 
#endif

		/* Test the new vector for a solution. */
		solutiontest(b[j-1], z, KMz, v1, v2, c1, nom, denom,flag);

		lllfp(b,mu,c,N,/*bs,*/j-1,h,z,delta,KMz,v1,v2,c1,nom,denom,flag);

/*
   If $N_{h-1} < 0.5$  the vectors in $b$ were linear dependent.
   Since at most one of the vectors $b_j,\ldots,b_{h}$ was linear dependent
   and is now zero at position $h$, swap the vectors back.
   Otherwise the vectors weren't linear dependent. This is the case if a
   multiple of $b_{s-1}$ was found in \|enumerate|.
*/
		if (fabs(N[h-1])<0.5) {
			swapvl = b[h-1];
			for (i=h;i<=s;i++) b[i-1] = b[i];
			b[s] = swapvl;
		}
		zaehler = 0;
	}
	else {
#if VERBOSE > 2
			printf("Not down %d  %d \n",j,k);  
#endif
		zaehler += 1;
/* Merely size-reduce the vectors. */
	}
	lllfp(b,mu,c,N,/*bs,*/h-1,h,z,0.0,KMz,v1,v2,c1,nom,denom,flag);
     }
     /* end \|while|  */
   }

/* Free the occupied memory. */
     lllfree(mu,c,N,/*bs,*/s);
     free(u);
     
     return(regulaerer_Durchlauf);
}

/* 
\section{Brute force enumeration}

\subsection{Final enumeration version of Schnorr, Euchner}

{\bf New} Enumerate in depth first search. Version without goto's 
as it is described in H.~H.~Hoerner's diploma thesis.

All important arrays are \|DOUBLE| arrays.

29.10.1996: \|tmax| is introduced.
But this version is not longer used. The next procedure is 
faster.
*/
DOUBLE enumerate_final(DOUBLE **mu, DOUBLE *c, long *u, int s, int j, int k)
{
	DOUBLE cd;
	DOUBLE *y;
	DOUBLE *cs;
	DOUBLE *us;
	DOUBLE dum;
	int t,i,tmax;
	int utest;
	
/* allocate memory for \|us| */
     us=(DOUBLE*)calloc(s+1,sizeof(DOUBLE));
     if (us == NULL) return(zu_wenig_Speicherplatz);

/* allocate memory for \|cs| */
     cs=(DOUBLE*)calloc(s+1,sizeof(DOUBLE));
     if (cs == NULL) return(zu_wenig_Speicherplatz);

/* allocate memory for \|y| */
     y=(DOUBLE*)calloc(s+1,sizeof(DOUBLE));
     if (y == NULL) return(zu_wenig_Speicherplatz);

/* Initiation */

     cd = fabs(c[j]);
     cs[k+1] = 0.0;

     /* Now we start from $t=j$ instead of $t=k$. */
     t = j;
     tmax = t;

     for (i=j+1;i<=k;i++) {
     		u[i] = 0;
     		us[i] = 0.0;
     		cs[i] = 0.0;
     		y[i] = 0.0;
     }

     u[j] = 1;
     y[t] = 0.0;
	
     /* *< us[t] = ceil(-SQRT(cd/c[t])); >* */
     /* With $t=j$ the last line is reduced to: */
     us[t] = -1.0;

/* The search loop. */     
     do {
     	     dum = us[t] + y[t];
	     cs[t] = cs[t+1] + dum*dum*c[t];
	     if (cs[t]<cd) {
		if (t>j) {
			/* Step forward. */
			t--;
		        for (i=t+1,dum=0.0;i<=tmax;i++) dum += us[i]*mu[i][t];
		        y[t] = dum;
		        us[t] = ceil(-y[t]-SQRT((cd - cs[t+1])/c[t]));
			continue;
		}
		else {
			/* Now we are at a leave of the enumeration tree. 
			   We have to test if we found a nontrivial 
			   solution.
			*/
			utest = 0;
			for (i=j;(i<=k);i++) if (us[i]!=0.0) { 
							utest = 1;
							break;
					     }
			if (utest==1) {
				cd = cs[j];
				for (i=j;i<=k;i++) u[i] = (long)us[i];			
			}
		}
	     }
	     else {
		/* Step back:
		
		   This \|if|-clause is new. 
		   \|us[t] + y[t]| is still stored in \|dum|.
		*/	     
	     	/* *< if (us[t] + y[t] >= -0.5) { >* */
	     	
	     	if (dum >= -0.5) {
	     		t++;
	     		if (tmax<t) tmax = t;
	     	}
	     }
	     if (t<=k) 	us[t] += 1;
     } while (t<=k);                                                                          
                                                  
/* Free the allocated memory. */
     free (us);
     free (cs);
     free (y);
     return (cd);
}

/* 
\subsection{Pruned Gauss-Enumeration}

{\bf New} Enumerate in depth first search. Pruned Version without goto's 
as it is described in H.~H.~Hoerner's diploma thesis.
*/
DOUBLE enumerate(DOUBLE **mu, DOUBLE *c, long *u, int s, int j, int k, int p)
{
	DOUBLE cd;
	DOUBLE *y;
	DOUBLE *cs;
	DOUBLE dum;
	DOUBLE *eta;

	long *us;
	int t,i;
	int utest;
	
	long *delta;
	long *d;
	long *v;
	int tmax;

/* Blocks with blocksize lower than \|SCHNITT| are completely enumerated,
   otherwise the enumeration is pruned. */
	int SCHNITT = 10;
	
/* $2^{-p}$ is the probability that lattice points are lost. */
	DOUBLE pi = 3.141592653589793238462643383;
	DOUBLE e = 2.718281828459045235360287471;
	
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

/* allocate memory for \|eta| */
     eta=(DOUBLE*)calloc(s+1,sizeof(DOUBLE));
     if (eta == NULL) return(zu_wenig_Speicherplatz);

/* allocate memory for \|v| */
     v=(long*)calloc(s+1,sizeof(long));
     if (v == NULL) return(zu_wenig_Speicherplatz);

/* Initiation */
     for (i=j;i<=k;i++) {
     	cs[i] = 0.0;
     	y[i] = 0.0;
     	u[i] = 0;
     	us[i] = 0;
     	v[i] = 0;
     	delta[i] = 0;
     	d[i] = 1;
     }
     	
     /* Now we start from $t=j$ instead of $t=k$. */
     t = j;
     tmax = t;
     cd = fabs(c[j]);
     us[j] = 1;
     u[j] = 1;
     cs[k+1] = 0.0;
     eta[j] = 0.0;

/* $\eta_t$ is precomputed, where
$$ \eta_t \approx {t\over 2\cdot \pi\cdot  e}
\biggl(\pi\cdot  t\cdot  p^2\cdot \prod_{s=1}^tc_s\biggr)^{1/t}.$$
In our case $s$ runs from $i$ to $k$ and $i=j+1,\ldots,k$. Therefore
$t=i-j$.
*/
     if (k-j < SCHNITT) {
     	for (i=j+1;i<=k;i++) eta[i] = 0.0;
     }
     else {
     	dum = c[j];
     	for (i=j+1;i<=k;i++) {
     		eta[i] = (i-j)*dum*exp( (log(pi*(i-j))-2.0*(DOUBLE)p*log(2.0)) /(i-j))/(2.0*pi*e);
#if VERBOSE > 2
 		printf("%13.10f\n",(double)eta[i]);
#endif
     		dum = exp( log(dum)*(i-j)/(i+1-j) + log(c[i])/(i+1-j));
     	}
     }

/* The search loop. */     
     do {
#if VERBOSE > 3
       	     printf("enumerate: Start at %d, us[%d]: %ld\n",j,t,us[t]);
#endif

     	     dum = us[t] + y[t];
	     cs[t] = cs[t+1] + dum*dum*c[t];

/* \|eta| controls the pruning. */
	     if (cs[t]<cd-eta[t]) {
		if (t>j) {
			t--;
			delta[t] = 0;
		        dum = 0.0;
		        for (i=t+1;i<=tmax;i++) dum += us[i]*mu[i][t];
		        y[t] = dum;
		        v[t] = (long)(ceil(-y[t] - 0.5));
		        us[t] = v[t];
			if (v[t] > -y[t]) d[t] = -1;
			else d[t] = 1;
		}
		else {
			utest = 0;
			for (i=j;i<=k;i++) if (us[i]!=0) { 
							utest = 1;
							break;
					     }
			if (utest==1) {
				cd = cs[j];
				for (i=j;i<=k;i++) u[i] = us[i];
			}
			/* We go straight on back because of possible rounding 
			   errors during the comparison \|cs[t]<cd-eta[t]|. 
			*/			
 			goto stepup;
		}
	     }
	     else {
stepup:		t++;
		if (tmax<t) tmax = t;
		if (t<tmax) delta[t] *= -1;
		if (delta[t]*d[t]>=0) delta[t] += d[t];
		us[t] = v[t] + delta[t];
	     }
     } while (t<=k);                                                                          
                                                  
/* Free the allocated memory. */
     free (us);
     free (cs);
     free (y);
     free (delta);
     free (d);
     free (eta);
     free (v);
     return (cd);
}

/* The maximum norm enumeration: 

\includefile{enumall.c}
*/
/* \endc */
