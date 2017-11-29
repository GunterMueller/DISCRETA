/*
% $Header: /usr/local/cvsroot/designs/designl.c,v 1.15 1997/12/03 17:40:23 alfred Exp $

\input cnoweb
\title{LLL-reduction on Kramer-Mesner-type matrices}
\job{designl.c}
\synopsis{Written by Alfred Wassermann, 10.1.95.

Sparse structure added 6.6.1995.

The program needs the following command line parameters:
\item{(1)} Scaling parameter $c_0$, e.g. $c_0=30$.
\item{(2)} Scaling parameter $c_1$, e.g. $c_1=10$. Now it
           has to be eqaul to $lambda$.
Old:           
\item{(3)} Scaling parameter $c_2$, e.g. $c_2=1$.
\item{(4)} $\beta$, length of the Korkine-Zolotarev blocks. e.g. $\beta = 40$.
Old:
\item{(5)} Number of loops (shuffling $\to$  bkz-reduction $\to$ test for
	   solution $\to$ Kreher-Radziszowski $\to$ shuffling), e.g. 50.
Old:
\item{(6)} Name of the log file, where intermediate output is written. 

Old: (17.9.1997)
22.2.1995: additional parameters:
$\gamma \cdot n$ is the estimated number of 1's in the solution vector, 
where $n$ is the number of lattice basis vectors.
\item{(7)} nominator of $\gamma$.
\item{(8)} denominator of $\gamma$.

13.6.1995: One more additional parameter:
\item{(9)} $p$, the pruning parameter. $p=10$ means fast pruning, $p=18$ means
    nearly complete enumeration.
    
24.1.97: It's possible to use an arbitrary right hand side:
The input file has to contain a 1 after the number of columns and of
course the column to the right of the matrix.
The number of columns is still the number of columns of the 
Kramer Mesner matrix, i.e. left hand side.

17.9.1997:
Selection of columns is enabled. If there is a flag n>0 on the command line,
the standard input is searched for a line
SELECTION n
and a line where the selected columns are indicated by 1.
No selection is done when there is no flag or flag$=0$.
		    
Normally $c_2$ is chosen to be 1. $30\leq \beta \leq 50$ is good choice. 
$c_0$ and $c_1$ should be chosen appropiately, so that the number which 
appears after `Analyse' decreases.
If nothing else is known $1/2$ would be a good choice for $\gamma$.

Starting the program on a unix system would look like:

\centerline{\|designl 30 10 1 40 50 designl.log 1 2 <KM.dat|}

where the file \|KM.dat| contains the Kramer-Mesner matrix.
It has to have the following form: The first lines of the  input file contain
some comments and have to be preceded by a \% in the first column.
The first line after the comments contains the number of rows and the number of
columns of the Kramer-Mesner matrix separated by a space.
Then the matrix follows, where the columns have to be separated by spaces.

{\bf There are no tests for wrong input data format!}

For example:

\|\%KM_PGGL_2_2h5_k5.txt|

\|\%This is a comment|

\|3 13|

\|4 8 4 4 4 4 0 0 0 0 0 0 0|

\|2 3 5 4 1 2 3 1 2 1 1 2 1|

\|0 5 0 5 5 0 0 0 5 1 1 5 1|


The algorithm searches for $\{0,1\}$ solutions of
Kramer-Mesner-type matrices via blockwise Korkine-Zolotarev reduction of
the following lattice (where it is estimated that the fraction of
1's in the solution vector is ${\gamma_n\over \gamma_d}$):
$$\pmatrix{ &    & &   c_0    & 0      & 0 \cr
	    & c_0KM & & \vdots & \vdots & \vdots \cr
	    &    & &  c_0     & 0      & 0 \cr
	  c_1\gamma_d &    & &  0     & c_1\gamma_n      & 0 \cr
	    & \ddots && \vdots & \vdots & \vdots \cr
	    &    &c_1\gamma_d&  0     & c_1\gamma_n      & 0 \cr
	  0  & \ldots   &0 &  c_2     & 0      & 0 \cr
	  0  & \ldots   &0 &  0     & 1      & 0 \cr}
$$

See Kreher and Radziszowski, ''Finding Simple $t$-Designs by Basis Reduction'',
{\it Congressus Numerantium} {\bf 55} (1986) 235--244.

\copyright\ \uppercase\expandafter{\romannumeral\year} by Alfred Wassermann.}

After reading the matrix with \|design_eingabe| and scaling
the resulting lattice with \|design_scale| we repeat the 
following steps:

\item{$\bullet$} The columns of the lattice are shuffled,
\item{$\bullet$} the columns which have a nonzero last row are written
	to the left of the lattice,
\item{$\bullet$} Blockwise Korkine Zolotarev reduction
\item{$\bullet$} Test if we found a solution
\item{$\bullet$} Simple enhancement of Kreher and Radziszowski.

The algorithm uses the subroutines in \|lllfpl.c| and
the include file \|llll.h|. 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "llll.h"

int selection;
int *permutation;
int inequalities;
int slackvar;
int c1_global;
int s_sel;
long stopafter;
long stoploops;
long initialguess;

long maximal_number_of_blocks;
long midpoint;
long grouporder;
long lambda;

long firstlinefactor = 10;

/* Some prototypes. */
void design_scale(COEFF **b, int s, int z, int KMz, long c0, long c1);
int design_eingabe(COEFF ***b, long **x, int *s, int *z, int *KMs, int *KMz, long nom, long denom, int selection, int **permutation, int *s_sel, char *filename);

long bin (long n, long k) {
	int i;
	long b;
	
	if ((n<k)||(n<0)||(k<0)) {
		printf("n=%ld and k=%ld not correct for bin\n", n,k);
		exit(0);
	}
	b=1;
	for(i=0;i<k;i++) {
		b *= (n-i);
		b /= (i+1);
	}
	return b;
}


/* \section{Main program} */
int main(int argc, char *argv[])
{
	int i,j;
	COEFF **b;
	long *x;
	int s,z,orgs;
	int KMs,KMz;

	long c0,c1;
	COEFF *swapvl;
	long v1,v2;
	long nom, denom;
	
	int flag;
	char *filename;

/* Here is the input of the command line parameters tested. We need exactly 9
   parameters. 
*/
	if ((argc<6)||(argc>7)) {
		printf("Wrong number of parameters in command line input!\n");
		printf("design c0 c1 beta p filename [SEL]\n");
		printf("design 1000 lambda 120 14 KM.dat 1\n");
		return 1;
	}

/* Test if they are well chosen. There a more sophisticated test will be 
   welcome. */
	for (i=1;i<=4;i++) {
		if ((atoi(argv[i])<=0)) {
			printf("Something is wrong with input parameter %d in the command line!\n",i);
			printf("design c0 c1 beta p filename [SEL]\n");
			printf("design 1000 lambda 120 14 KM.dat 1\n");
			return 1;
		}
	}
	
/*
The constants $c_0$, $c_1$, $c_2$. These constants are controlling
very sensitive the behaviour of the \|bkz|-algorithm.
*/
/*	c0 = atoi(argv[1]); */
	c0 = 1;

	initialguess = atoi(argv[1]);
	c1 = atoi(argv[2]);
	lambda = c1;

	nom = 1;
	denom = 2;
	filename = argv[5];
	if (argc>6) {
		selection = atoi(argv[6]);
		if (selection<0) selection = 0;
	} else {
		selection = 0;
	}
	s_sel = 0;

/*   Read $b$ and $x$. $b$ has $s-1$ columns plus one zero column
     for LLL reduction.
*/
	nom   = 1 * (initialguess/2+1);
	denom = 2 * (initialguess/2+1);

	flag = design_eingabe(&b,&x,&s,&z,&KMs,&KMz,nom,denom,selection,&permutation,&s_sel,filename);
	if (flag==1) c1 = 1;
	if (inequalities==1) c1 = c1_global;

	orgs = s;

/* The $(0,\ldots,0,1,\ldots,1)^t$ column which we want to
approximate is written to position $0$.
*/
   swapvl = b[s-2];
	for (i=s-2;i>0;i--) b[i] = b[i-1];
	b[0] = swapvl;

	/*   Scale the matrix $b$ */
/*	design_scale(b,s,z,KMz,c0,c1);*/

	
	for (i=0;i<s;i++) b[i][1].c *= firstlinefactor;

/* We will have found a solution if a column consists of zeros in
the first rows (the Kramer-Mesner part of the lattice), $\pm 1$ in the
last row $(z-1)$ and only either $v_1$ or $v_2$ in all other rows except
row $z-2$. If row $z-2$ is different from zero we found a
non trivial solution.
*/
	v1 = c1;
	v2 = -c1;

	/*   Basis reduction 

	First we reduce fast, to get rid of some columns. */
#if VERBOSE > 2
	printf("Start of first reduction\n");
#endif
	lll(b,s-1,z,KMz,v1,v2,c1,nom,denom,0);
#if VERBOSE > 2
	printf("end of first reduction\n");
	ausgabe(b,s,z); 
#endif

	/* Nevertheless we test the basis for solutions */
	for (j=0;j<s-1;j++) 
		solutiontest(b[j], z, KMz, v1, v2, c1, nom, denom,0);

	/* Delete the unnecessary columns and rows. 
	   This was done by writing in a temporary file lll.output.
	   Now it is done directly.
	*/
	/* *< 
		putlattice(b,s,z,KMz,KMs); 
		getlattice(b,&s,&z,1); 
	>* */

/*	cutlattice(b,&s,&z,KMz); */

	/* Reduce with blockwise Korkine Zolotarev reduction. */  
#if VERBOSE > 2
	printf("Start of second reduction\n");
#endif
	bkz(b,s-1,z,LLLCONST_HIGH,atoi(argv[3]),KMz,v1,v2,c1,nom,denom,0,atoi(argv[4])); 
#if VERBOSE > 1
	printf("end of second reduction\n");
	ausgabe(b,s,z); 
#endif
	for (i=0;i<s;i++) b[i][1].c /= firstlinefactor;

	for (i=0;i<s;i++) printf("%ld ",b[i][1].c);
	printf("\n\n");

	/* Complete enumeration of all solutions. */  
	enuminf(b,s-1,z,KMs,KMz,v1,v2,/*c0,c1,*/nom,denom,selection,permutation); 

	/* Free the allocated memory */
	for (i=0;i<s;i++) free(b[i]);
	free(b);
	free(x);
	if (selection>0) free(permutation);
	printf("finished.\n"); fflush(stdout);

	return 0;
	/* finito!!! */
}

/* \section{Special functions}
\subsection{Scaling function}
The Kramer-Mesner matrix is multiplied by $c_0$.
The identity matrix which is appended to the Kramer-Mesner matrix is
scaled with $c_1$. The last two rows are multiplied by $c_2$.
*/
void design_scale(COEFF **b, int s, int z, int KMz, long c0, long c1)
{
	int i,j;
	for (i=0;i<s;i++) {
		for (j=1;j<=KMz;j++) b[i][j].c *= c0;
		if (inequalities==1) {
			if (i==0) {
				for (j=KMz+1;j<z-1;j++) b[i][j].c *= c1;
			} else {
				for (j=KMz+1;j<z-1-slackvar;j++) b[i][j].c *= c1;
			}
		} else {
			for (j=KMz+1;j<z-1;j++) b[i][j].c *= c1;
		}
		b[i][z].c *= c1;
	}
	return;
}
/* 
\subsection{Input function}
Input of the Kramer-Mesner-matrix.
The first lines of the input file can contain some comments
which must be preceded by a \% in the first column.
After the comments the following line contains the number of rows and 
the number of columns of the Kramer-Mesner matrix.
There is nearly no test for wrong data format.

Append the $I$-matrix,
the $\lambda$-column, $(0,\ldots,0,1\ldots,1)^t$ and a $0$-column.
The matrix $I$ is multiplied by the denominator of $\gamma$, say $\gamma_d$.
The appended vector is multiplied by the nominator $\gamma_n$.
The lattice has the following form:
$$\pmatrix{ &    & &   1    & 0      & 0 \cr
	    & KM & & \vdots & \vdots & \vdots \cr
	    &    & &  1     & 0      & 0 \cr
	  \gamma_d &    & &  0     & \gamma_n      & 0 \cr
	    & \ddots && \vdots & \vdots & \vdots \cr
	    &    &\gamma_d&  0     & \gamma_n      & 0 \cr
	  0  & \ldots   &0 &  1     & 0      & 0 \cr
	  0  & \ldots   &0 &  0     & 1      & 0 \cr}
$$
The last column is needed in \|bkz|.

Returns 1, if there is an arbitrary right hand side. 
0 otherwise.
*/

int design_eingabe(COEFF ***b,long **x, int *s, int *z, int *KMs, int *KMz, long nom, long denom, int selection, int **permutation, int *s_sel, char *filename)
{
#define zlength 10000

	int i,i0,j0,j,l,flag;
	char zeile[zlength];
	char select[100];
	COEFF *swap;
	int swap2;
	char *zp;
	FILE *txt;
	int KMs_org;
	
	int t,v,k;
	long w;

   txt = fopen(filename, "r");
	sprintf(select,"%% with RHS");
	inequalities = 0;
	stopafter = 0;
	stoploops = 0;
	do {
		fgets(zeile,zlength,txt);  
		if (strstr(zeile,"% degree:")!=NULL) {
			sscanf(zeile,"%% degree: %d",&v);
		}
		if (strstr(zeile,"% t, k:")!=NULL) {
			fgets(zeile,zlength,txt);
			sscanf(zeile,"%% %d %d",&t,&k);
		}
      if (strstr(zeile,"% order:")!=NULL) {
			sscanf(zeile,"%% order: %ld",&grouporder);
		}
		
      if (strstr(zeile,select)!=NULL) inequalities = 1;
      if (strstr(zeile,"% stopafter")!=NULL) {
			sscanf(zeile,"%% stopafter %ld",&stopafter);
		}
      if (strstr(zeile,"% stoploops")!=NULL) {
			sscanf(zeile,"%% stoploops %ld",&stoploops);
		}
	}
	while (zeile[0]=='%');

	maximal_number_of_blocks =  lambda * bin(v,t) / bin(k,t);
/*	maximal_number_of_blocks =  84; */
	midpoint = maximal_number_of_blocks-nom+1;

	printf("%d-(%d,%d,%ld) designs\n",t,v,k,lambda); 
	printf("Maximal number of blocks: %ld\n",maximal_number_of_blocks); 
	printf("Midpoint: %ld, nom: %ld\n",midpoint, nom); 
	fflush(stdout);
	
   /* \|flag|$=1$ means arbitrary right hand side. */
	flag = 0;       
	sscanf(zeile,"%d%d%d",&(*KMz),&(*KMs),&flag); 
	if (flag!=1) sscanf(zeile,"%d%d",&(*KMz),&(*KMs));

	KMs_org = (*KMs);
	if (inequalities==1) {
		(*KMs) += (*KMz);
	}
/*------------------------------------*/	
/* +1 for packings */	
/*	(*z) = (*KMz)+(*KMs)+2 +1;*/
/*------------------------------------*/	
	(*z) = (*KMz)+(*KMs)+2 ;
	/* Two additional columns: the $\lambda$-column and a zero-column. */
	(*s) = (*KMs)+2;


/* allocation */
	(*x) = (long*)calloc((*s),sizeof(long));
	for (i=0;i<(*s);i++) (*x)[i] = 0;

	(*b) = (COEFF**)calloc((*s),sizeof(long*));
	for(j=0;j<(*s);j++) {
		(*b)[j] = (COEFF*)calloc((*z)+1,sizeof(COEFF));
		for (i=0;i<=(*z);i++) (*b)[j][i].c = 0; 
	}
	
/* Read the Kramer-Mesner-matrix */
	for (j=0;j<(*KMz);j++) {
		for (i=0;i<KMs_org;i++) {
			fscanf(txt,"%ld",&((*b)[i][j+1+1].c));
			(*b)[i][j+1+1].c *= denom;
		}
		(*b)[(*KMs)][j+1+1].c = lambda*nom; 
		(*b)[(*KMs)+1][j+1+1].c = 0;            
	}

/* 
   Append the other columns and lines.
	Meanwhile, \|denom|=2, \|nom|=1.
*/

	for (j=(*KMz)+1;j<(*z);j++) {
		for (i=0;i<(*s)-1;i++) (*b)[i][j+1].c = 0;
		(*b)[j-(*KMz)-1][j+1].c = denom;
		(*b)[(*s)-2][j+1].c = nom;
		(*b)[(*s)-1][j+1].c = 0;
	}
	(*b)[(*s)-2][(*z)-1].c = 1*nom; 
/*	(*b)[(*s)-2][(*z)-1].c = 0; */
	(*b)[(*s)-2][(*z)].c = 1*nom;
/* 
	Read the weights
*/
	sprintf(select,"STABILIZER-ORDER-K-SETS");
	do {
		zp=fgets(zeile,zlength,txt);  
	} while ((zp!=NULL)&&(strstr(zeile,select)==NULL));

	if (zp==NULL) {
		printf("%s not found\n",select); fflush(stdout);
		fclose(txt); 
		exit(0);
	}
	for (i=0;i<(*KMs);i++) {
		fscanf(txt,"%ld",&w);
		(*b)[i][0+1].c = grouporder/w;
/*		(*b)[i][0+1].c = 0; */
	}
	(*b)[(*KMs)][0+1].c = midpoint; 

/* 
	Handle the selection vector
*/
	(*s_sel) = (*s);
	if (selection>0) {
		sprintf(select,"SELECTION %d",selection);
		do {
			zp=fgets(zeile,zlength,txt);  
		}
		while ((zp!=NULL)&&(strstr(zeile,select)==NULL));

		if (zp==NULL) {
			printf("%s not found\n",select);
			fclose(txt); 
			exit(0);
		}
		zp=fgets(zeile,zlength,txt); 
		fclose(txt); 
		if (zp==NULL) {
			printf("Selection vector not found\n");
			exit(0);
		}
		printf("%s:\n",select);
		printf("%s\n",zeile);
		
		(*permutation) = (int*)calloc((*s),sizeof(int));
		i0=0;
		zp = zeile;
		for (i=0;i<(*s_sel);i++) (*permutation)[i]=i;
		for (i=0;i<(*s_sel)-3;i++) {
			sscanf(zp,"%d",&l);
			zp=&(zp[2]);
			/* columns corresponding to 0 are deleted */
			if (l==0) {
				swap = (*b)[i0];
				swap2 = (*permutation)[i0];
				for (j=i0+1;j<(*s_sel);j++) {
					(*b)[j-1] = (*b)[j];
					(*permutation)[j-1] = (*permutation)[j];
				}
				(*b)[(*s_sel)-1] = swap;
				(*permutation)[(*s_sel)-1] = swap2;
				(*s)--;
				for (j0=0;j0<(*s_sel);j0++) {
					for (j=(*KMz)+i0+1;j<(*z);j++) {
						(*b)[j0][j].c = (*b)[j0][j+1].c;
					}
				}
				(*z)--;
			} else {
				i0++;
			}
		}
		(*KMs) = (*s)-3;
	}

/* Find the nonzero entries. */
	for (i=0;i<(*s_sel)-1;i++) coeffinit((*b)[i],(*z));
	
	return flag;
}

/* 
\subsection{Test for solutions}
This test can be executed whenever a lattice vector is modified.
*/

int solutiontest(COEFF *b, int z, int KMz, long v1, long v2, long c1, long nom, long denom, int flag)
{
	int i;
	int plus;
	int F;
	int upper;
	int KMzlocal;
	long noblocks;
	
	long u;
	int i0,l;
	
/* Test if we found a solution. But we are only interested in nontrivial 
   designs ($\lambda\neq 0$).
*/
	F = 0;
		
/* \|labs(b[j][z-1])==1.0| means the solution is inhomogenous.
   \|b[j][z-2]| indicates a nontrivial solution and \|plus==0| indicates
   a solution.
*/   
	if ((labs(b[z].c)!=nom) || (b[z-1].c==0)) return 0;
	plus = 0;

/* We are still working with the complete lattice. The Kramer Mesner part
    is not deleted.
*/      
/*
	if (flag==0) for(i=1;i<=KMzlocal;i++) plus += b[i].c * b[i].c;
*/
	
/* if \|flag|${} = 1$ \|KMz| is locally changed. */     
	KMzlocal = 1;
		
	if (plus==0) {
		upper = z;
		F = 1;
		for (i=KMzlocal+1;i<KMz+2;i++) {
			if (labs(b[i].c)>nom) {
				F = 0;
				break;
			}
		}
		if (F) {
			for (i=KMz+2;i<upper;i++) {
				if (labs(b[i].c)!=nom) {
					F = 0;
					break;
				}
			}
		}
			
		/* We found a solution, Yipeeeh! */
		if (F==1) {
			printf(" ");
			if (!selection) {		
/*
				for (i=1;i<upper;i++) {
					printf("%ld ",b[i].c);
				}
				printf("\n");
*/				
				for (i=KMz+1+1;i<upper;i++) {
					u = labs(b[i].c-b[z].c)/denom/c1;
					printf("%ld ",u);
				}
			} else {

/* SELECTION and RHS are not yet synchronized */
				i0 = 0;
				for (i=KMzlocal;i<KMzlocal+s_sel-3;i++) {
					u = 0;
					for (l=i0;l<=i-KMzlocal;l++) {
						if (permutation[l]==i-KMzlocal) {
							u = labs(b[i0+1].c-b[z].c*nom)/denom/c1;
							i0++;
							break;
						}
					}
					printf("%ld ",u);
				}
			}
			if (b[1].c>0) {
				noblocks = midpoint-b[1].c;
			} else {
				noblocks = midpoint+b[1].c;
			}
			printf(" blocks = %ld of %ld\n",noblocks,maximal_number_of_blocks);
			fflush(stdout);
		}
	}
	return F;
}

/* \endc */
