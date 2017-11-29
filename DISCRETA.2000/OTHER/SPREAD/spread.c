/*
  $Header: /usr/local/cvsroot/spread/spread.c,v 1.14 1998/03/30 16:59:04 alfred Exp $
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define zlength 16000
#define MAXINC 10000
#define VERBOSE 0

FILE *infile;
FILE *outfile;	
FILE *solfile;
int lambda;
int numeqn, numvar;
int maxlines,maxdepth;

int *Pi,*Li,*Si,**Lp,*Sel;
int *Partial;
/* int ***Mr, ***Mc; */
int **Mr, **Mc;
int *Mcnum, *Mrnum;
int *Pihist,*Lihist,*Sihist;
int *Pihistp,*Lihistp,*Sihistp;

long nosolutions;
long loops;
long stopafter;
long stoploops;

int selection;

void changeone(int level, int *list, int *hist, int *pointer, int pos, int sign) {
	if (sign==-1) {
		list[pos]--;
	} else {
		list[pos]++;
	}
	hist[pointer[level+1]] = pos;
	pointer[level+1]++;	
	return;
}

void changePi(int level, int *list, int *hist, int *pointer, int pos, int entry) {
	list[pos] -= entry;
	hist[pointer[level+1]] = pos;
	pointer[level+1]++;	
	hist[pointer[level+1]] = entry;
	pointer[level+1]++;	
	
	Partial[level] -= entry;
	
	return;
}

void restoreone(int level, int *list, int *hist, int *pointer, int sign) {
	int i;
	for (i=pointer[level];i<pointer[level+1];i++) {
		if (sign==-1) {
			list[hist[i]]++;
		} else {
			list[hist[i]]--;
		}
	}
	pointer[level+1]=pointer[level];
	return;
}

void restorePi(int level, int *list, int *hist, int *pointer) {
	int i;
	for (i=pointer[level];i<pointer[level+1];i+=2) {
			list[hist[i]] += hist[i+1];
			Partial[level] += hist[i+1];
	}
	pointer[level+1]=pointer[level];
	return;
}
/*-----------------------------------------------------------------------*/
/*
void search_point(int *numlines, int *p) {
	static int i, j, b;
	
	(*numlines) = MAXINC;
	(*p)=-1;
	for (i=0;i<numeqn;i++) if (Pi[i]>0) { 
		b = 0;
		for (j=Mrnum[i]-1;j>=0;j--) if (Li[Mr[i][j][1]]==1) {
			b++;
			if (b>=(*numlines)) break;
		}
		if (b<(*numlines)) {
			(*numlines) = b;
			(*p) = i;		
			if ((*numlines)<=1) break;	
		}	
	}
	return;
}	
void generate_Lp(int level, int p) {
	static int i,l;
	i = numlines-1;
	for (l=Mrnum[p]-1;l>=0;l--) if (Li[Mr[p][l][1]]==1) {
		Lp[level][i] = l;
		if (i==0) break;
		i--;
	}
	return;
	i = 0;
	for (l=0;l<Mrnum[p];l++) if (Li[Mr[p][l][1]]==1) {
		Lp[level][i] = l;
		i++;
	}
}	
*/
void spread(int level) {
	int i,j,p,l,k,li;
	int b;
	int numlines;
	int line;
	int point;
	
#if VERBOSE > 0
	printf("Level: %d\n",level);
#endif 
#if VERBOSE > 1
	for (i=0;i<numvar;i++) printf("%d ",Li[i]); printf("\n");
	for (i=0;i<numvar;i++) printf("%d ",Si[i]); printf("\n");
	for (i=0;i<numeqn;i++) printf("%d ",Pi[i]); printf("\n");
#endif 
	if (level==0) {
		p = lambda*numeqn;
		Partial[level] = p;
	} else {
		Partial[level] = Partial[level-1];
		p = Partial[level];
	}

#if VERBOSE > 0
if (loops%20000==0) {
	printf("Loops: %ld, level: %d\n",loops,level);
	j = 0;
	for (i=0;i<numvar;i++) {
		if (Si[i]==1) {
			printf("1");
		} else {
			if (Li[i]==0) {
				printf("0");
			} else {
				printf(".");
				j++;
			}
		}
	}
	printf("\nDots: %d, level: %d, Points: %d\n\n",j,level,Partial[level]);
	fflush(stdout);
}	
#endif

/* If Pi is empty, print Si as solution */
	if (p<=0) {
#if defined(SHORTOUTPUT)		
/*		for (i=0;i<numvar;i++) if (Si[i]==1) printf("%d ",i); printf("\n"); fflush(stdout); */
		for (i=0;i<numvar;i++) if (Si[i]==1) fprintf(solfile,"%d ",i); fprintf(solfile,"\n");
#else
		for (i=0;i<numvar;i++) printf("%d",Si[i]); printf("\n"); fflush(stdout);
		for (i=0;i<numvar;i++) fprintf(solfile,"%d",Si[i]); fprintf(solfile,"\n");		
#endif		
		fflush(solfile);
		nosolutions++;
		if ((stopafter>0)&&(stopafter<=nosolutions)) {
			printf("Loops: %ld\n",loops);
			printf("Stop after %ld solutions.\n",nosolutions);
			exit(1);
		}
		return;		 
	}

	loops++;
	if ((stoploops>0)&&(stoploops<loops)) {
		printf("Stop after %ld loops\n",loops);
		printf("Number of solutions: %ld\n",nosolutions);
		exit(1);
	}
/* 
	There are still points in Pi, we have to go on.
   Find a point p in Pi incident with the fewest lines Li(p) in Li. 
*/
/*
	search_point(&numlines, &p);
*/	

	numlines = MAXINC;
	for (i=0;i<numeqn;i++) if (Pi[i]>0) { 
		b = 0;
		for (j=Mrnum[i]-1;j>=0;j--) if (Li[Mr[i][j]]==1) {
			b++;
			if (b>=numlines) break;
		}
		if (b<numlines) {
			numlines = b;
			p = i;  	 
			if (numlines<=0) break;	 
		}	
	}

	
	if (numlines==0) {
#if VERBOSE > 2
		printf("No line to choose -> return\n");
#endif
		return;
	}		


/* Generate L(p) with numlines entries */
/*
	generate_Lp(level, p);
*/		
	i = numlines-1;
	for (l=Mrnum[p]-1;l>=0;l--) if (Li[Mr[p][l]]==1) {
		Lp[level][i] = l;
		if (i==0) break;
		i--;
	}

#if VERBOSE > 1
	printf("Lines in L(%d): ",p);
	for (l=0;l<numlines;l++) printf("%d-%d ",Lp[level][l],Mr[p][l]);
	printf("\n");
#endif	

/* For each l in Lp do */	
	for (l=0;l<numlines;l++) {
		Pihistp[level+1] = Pihistp[level];
		Lihistp[level+1] = Lihistp[level];
		Sihistp[level+1] = Sihistp[level];

		line = Mr[p][Lp[level][l]];
#if VERBOSE > 1
		printf("Choose line: %d, Number of points: %d\n",line,Mcnum[line]);
		fflush(stdout);
#endif
		/* add line to the partial spread */
		changeone(level,Si,Sihist,Sihistp,line,1); 
		
		if (lambda==1) {
		/* lambda = 1 version */
			for (i=Mcnum[line]-1;i>=0;i--) {
				point = Mc[line][i<<1];
				if (Pi[point]>0) {
					/* delete point */
					changePi(level,Pi,Pihist,Pihistp,point,1); 
					/* delete lines incident to this point */
					for (j=Mrnum[point]-1;j>=0;j--) if (Li[Mr[point][j]]==1) {
						changeone(level,Li,Lihist,Lihistp,Mr[point][j],-1); 
					}	
				}
			}
		} else {
		/* lambda > 1 version */
			/* delete all previous chosen lines to avoid redundancy */
			for (i=l;i>=0;i--) {
				k = Lp[level][i];
				if (Li[Mr[p][k]]==1) changeone(level,Li,Lihist,Lihistp,Mr[p][k],-1);
			}			
			/* subtract points */			
			for (i=Mcnum[line]-1;i>=0;i--) {
				point = Mc[line][i<<1];
				changePi(level,Pi,Pihist,Pihistp,point,Mc[line][(i<<1)+1]); 
			}
			/* Search for lines incident with the above point which are no
			   longer possible.
			 */
			for (i=Mcnum[line]-1;i>=0;i--) {
				point = Mc[line][i<<1];
				for (j=Mrnum[point]-1;j>=0;j--) {
					li = Mr[point][j];
					if (Li[li]==1) {
						for (k=Mcnum[li]-1;k>=0;k--) if (Mc[li][(k<<1)+1]>Pi[Mc[li][k<<1]]) {
							changeone(level,Li,Lihist,Lihistp,li,-1); 
							break;
						}
					}
				}	
			}
		}
		/* Recursion */
 		spread(level+1);
		
		/* Restore the old values before the recursion */
		restorePi(level,Pi,Pihist,Pihistp);
		restoreone(level,Li,Lihist,Lihistp,-1);
		restoreone(level,Si,Sihist,Sihistp,1);
	}
	return;
}
/*-----------------------------------------------------------------------*/
void main(int argc, char *argv[]) {
	int i,j,k,s,b,i_old;
	int b1;
	int j0;
	int level;
	char zeile[zlength];
	char *filename;
	char select[100];
	char *zp;
	
	if ((argc!=3)&&(argc!=4)) {
		printf("Wrong parameters! Correct are: lambda filename [SEL]\n");
		exit(0) ;
	}
 	lambda = atoi(argv[1]);
	filename = argv[2];
	nosolutions = 0;
	loops = 0;
	if (argc>3) {
		selection = atoi(argv[3]);
		if (selection<0) selection = 0;
	} else {
		selection = 0;
	}

/* Read the matrix and write the nonzero entries to 
   spreadtmp.txt. All columns with entries > lambda are deleted.
	Also the number of entries for each line and each column
	are counted in oreder to save space.
*/
	stopafter = -1;
	stoploops = -1;
	infile = fopen(filename,"r");
	do {
		fgets(zeile,zlength,infile);
      if (strstr(zeile,"% stopafter")!=NULL) {
			sscanf(zeile,"%% stopafter %ld",&stopafter);
		}
      if (strstr(zeile,"% stoploops")!=NULL) {
			sscanf(zeile,"%% stoploops %ld",&stoploops);
		}
	} while (zeile[0]=='%');
	sscanf(zeile,"%d%d",&numeqn,&numvar);

	printf("numvar: %d, numeqn: %d\n",numvar,numeqn);
	maxdepth = 2 + (numeqn*lambda<numvar)?numeqn*lambda:numvar;
		
	Pi = (int*)calloc(numeqn,sizeof(int));
	Si = (int*)calloc(numvar,sizeof(int));
	Li = (int*)calloc(numvar,sizeof(int));
	if (selection>0) {
		Sel = (int*)calloc(numvar,sizeof(int));
	}
	Pihist = (int*)calloc(2*maxdepth,sizeof(int));
	Sihist = (int*)calloc(numvar,sizeof(int));
	Lihist = (int*)calloc(numvar,sizeof(int));
	Sihistp = (int*)calloc(maxdepth,sizeof(int));
	Pihistp = (int*)calloc(maxdepth,sizeof(int));
	Lihistp = (int*)calloc(maxdepth,sizeof(int));

	Partial = (int*)calloc(maxdepth,sizeof(int));
	
	Mcnum = calloc(numvar,sizeof(int));
	Mrnum = calloc(numeqn,sizeof(int));
	
	for (i=0;i<numvar;i++) Si[i]=1;
	for (i=0;i<numvar;i++) Mcnum[i] = 0;
	for (i=0;i<numeqn;i++) Mrnum[i] = 0;
	maxlines = 0;
	outfile = fopen("spreadtmp.txt","w");
	for (i=0;i<numeqn;i++) {
		for (j=0;j<numvar;j++) {
			fscanf(infile,"%d",&s);
			if (s>0) {
				fprintf(outfile,"%d %d %d\n",i,j,s);
				Mcnum[j]++;
				Mrnum[i]++;
				if (s>lambda) Si[j] = 0;  
			}
		}
		if (Mrnum[i]>maxlines) maxlines = Mrnum[i];
	}
#if VERBOSE > 0
	printf("columns to keep:\n");
	for (i=0;i<numvar;i++) printf("%d ",Si[i]); printf("\n");
#endif	
	if (selection>0) {
		sprintf(select,"SELECTION %d",selection);
		do {
			zp=fgets(zeile,zlength,infile);  
		}
		while ((zp!=NULL)&&(strstr(zeile,select)==NULL));

		if (zp==NULL) {
			printf("%s not found\n",select);
			fclose(infile); 
			exit(0);
		}
		zp=fgets(zeile,zlength,infile); 
		if (zp==NULL) {
			printf("SELECTION vector not found\n");
			exit(0);
		}
		printf("%s:\n",select);
		printf("%s\n",zeile);
		
		zp = zeile;
		for (i=0;i<numvar;i++) {
			sscanf(zp,"%d",&k);
			zp=&(zp[2]);
			if (k==1) {
				Sel[i] = 1;
			} else {
				Sel[i] = 0;
			}
		}
	}
	fclose(infile);
	fclose(outfile);
/* 
	Now, the nonzero entries have been written into the file
   spreadtmp.txt. We allocate memory and read them again.
*/	

	b1 = 0;
	for (i=0;i<numeqn;i++) {
		b1 += Mrnum[i];
	}
#if VERBOSE > 0	
	printf("number of entries in column:\n");
	for (i=0;i<numvar;i++) printf("%d ",Mcnum[i]); printf("\n");
	printf("number of entries in row:\n");
	for (i=0;i<numeqn;i++) {
		printf("%d ",Mrnum[i]);
	}
	printf("\n");
	printf("Entries:: %d\n",b1);	
#endif	
	
	Lp = calloc(maxdepth,sizeof(int*));		
	for (i=0;i<maxdepth;i++) Lp[i] = calloc(maxlines,sizeof(int));

/*	
	Mr = calloc(numeqn,sizeof(int**));
	for (i=0;i<numeqn;i++) {
		Mr[i] = calloc(Mrnum[i],sizeof(int*));
		for (j=0;j<Mrnum[i];j++) Mr[i][j] = calloc(2,sizeof(int));
	}
	Mc = calloc(numvar,sizeof(int**));
	for (i=0;i<numvar;i++) {
		Mc[i] = calloc(Mcnum[i],sizeof(int*));
		for (j=0;j<Mcnum[i];j++) Mc[i][j] = calloc(2,sizeof(int));
	}
*/	
	Mr = calloc(numeqn,sizeof(int*));
	for (i=0;i<numeqn;i++) {
		Mr[i] = calloc(Mrnum[i],sizeof(int));
	}

	Mc = calloc(numvar,sizeof(int*));
	for (i=0;i<numvar;i++) {
		Mc[i] = calloc(2*Mcnum[i],sizeof(int));
	}

	infile = fopen("spreadtmp.txt","r");
	b=0;
	for (i=0;i<numeqn;i++) {
		b += Mrnum[i];
	}
#if VERBOSE > 0
	printf("number of entries in row again:\n");
	for (i=0;i<numeqn;i++) {
		printf("%d ",Mrnum[i]);
	}
	printf("\n");
	printf("Entries: %d\n",b);	
#endif

	for (i=0;i<numvar;i++) Li[i] = 0;
	j0 = 0;
	for (k=0;k<b;k++) {
		i_old = i;
		fscanf(infile,"%d%d%d",&i,&j,&s);
		if (i!=i_old) j0 = 0;
		if (Si[j]==1) {		
/*			Mr[i][(j0<<1)+1] = s; */
			Mr[i][j0] = j;
			j0++;
		
			Mc[j][(Li[j]<<1)+1] = s;
			Mc[j][Li[j]<<1] = i;
			Li[j]++;
		} else {
			Mcnum[j] = 0;
			Mrnum[i]--;
			printf("delete line %d\n",j);
			fflush(stdout);
		}
	}
	fclose(infile);

	printf("Input finished\n");
	fflush(stdout);
/* Now the matrix has been read again and it is available twice:
   via column number and via row number in Mc resp. Mr.

	Initially, the arrays for points and lines are cleared. 
	The partial spread is empty. 
*/	
	level = 0;
	for (i=0;i<numeqn;i++) Pi[i] = lambda;
	for (i=0;i<numvar;i++) if (Mcnum[i]>0) {
		Li[i] = 1;	
	} else {
		Li[i] = 0;	
	}
	for (i=0;i<numvar;i++) Si[i] = 0;
	if (selection>0) {
		for (i=0;i<numvar;i++) if (Sel[i]>0) {
			Li[i] = 1;	
		} else {
			Li[i] = 0;	
		}
	}
	k = 0;
	for (i=0;i<numvar;i++) if (Li[i]==1) k++;
	printf("Lines: %d, points: %d\n",k,numeqn); fflush(stdout);
	Pihistp[0] = 0;
	Lihistp[0] = 0;
	Sihistp[0] = 0;
	Pihistp[1] = 0;
	Lihistp[1] = 0;
	Sihistp[1] = 0;
#if VERBOSE > 0	
	printf("Start of recursion\n"); fflush(stdout);
#endif	
	solfile = fopen("solutions","w");

	spread(level);
	
	fclose(solfile);
	printf("Loops: %ld\n",loops);
	printf("Number of solutions: %ld\n",nosolutions);
	printf("finished.\n"); fflush(stdout);
	
/*
	printf("average number of fewest entries: %f\n",1.0*average/(double)averagen);
*/	
}
