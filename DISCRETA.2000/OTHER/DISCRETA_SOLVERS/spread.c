/*
  $Header: /usr/local/cvsroot/spread/spread.c,v 1.29 1999/10/04 11:55:58 alfred Exp $
*/
/******************************************************************
*
* Solve linear equations over the natural numbers 
* with 0/1 vectors.
*
* Author: Alfred Wassermann
* Date: 1998
*
* Literature:
*  R. Mathon,
*  "Searching for spreads and packings,"
*  Geometry, Combinatorial Designs and Related Structures
*  London Math. Society LNS 245
*******************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define FIX_FULL_ORBITS 0
/*#undef FIX_FULL_ORBITS*/

/* for timings: */
#include <sys/times.h>
#include <unistd.h>
#include <limits.h>
/* end timings */

/******************************
*
*  Flag for Random version:
*
*******************************/
#define RANDOM 0
#define STARTBACKTRACK 20
/*******************************/

#define zlength 31000
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
int **Mr, **Mc;
int *Mcnum, *Mrnum;
int *Pihist,*Lihist,*Sihist;
int *Pihistp,*Lihistp,*Sihistp;

long nosolutions;
long loops;
long stopafter;
long stoploops;

int selection;
int bestsofar;

#if RANDOM
int jump_back = 0;
#endif

#if defined(FIX_FULL_ORBITS)
int fullorbits;
long numberofblocks;
long order;
long *blocks;
int *Lifull;
int *Lishort;
long switchcapacity,blockstodo;
long numberswitch;
int numberfulls;
#endif
int switched=0;

int SILENT;

/******************************************************************
*  Functions to measure the runtime
*******************************************************************/
int user_time, time_0, time_1;
char timestring[256];

int os_ticks()
{
	struct tms tms_buffer;

	if (-1 == times(&tms_buffer))
		return(-1);
	return(tms_buffer.tms_utime);
}

int os_ticks_per_second()
{
	int clk_tck = 1;
	
	clk_tck = sysconf(_SC_CLK_TCK);
	/* printf("clk_tck = %ld\n", clk_tck); */
	return(clk_tck);
}

int os_ticks_to_hms_tps(int ticks, int tps, int *h, int *m, int *s)
{
	int l1;

	l1 = ticks / tps;
	*s = l1 % 60;
	l1 -= *s;
	l1 /= 60;
	*m = l1 % 60;
	l1 -= *m;
	l1 /= 60;
	*h = l1;
	return(1);
}

int os_ticks_to_hms(int ticks, int *h, int *m, int *s)
{
	os_ticks_to_hms_tps(ticks, os_ticks_per_second(), h, m, s);
	return(1);
}

void print_delta_time_tps(int l, int tps, char *str)
{
	int h, m, s;

	os_ticks_to_hms_tps(l, tps, &h, &m, &s);
	sprintf(str, "%d:%02d:%02d", h, m, s);
}

void print_delta_time(int l, char *str)
{
	print_delta_time_tps(l, os_ticks_per_second(), str);
}

/******************************************************************
*  Functions to update the lists.
*  The changes are stored in history-arrays and
*  can be restored fast.
*******************************************************************/
void change_one_in_list(int level1, int *list, int *hist, int *pointer, int pos, int sign) {
	if (sign==-1) {
		list[pos]--;
	} else {
		list[pos]++;
	}
	hist[pointer[level1]] = pos;
	pointer[level1]++;	
	return;
}

void change_Li(int level1, int *list, int *hist, int *pointer, int pos) {
	list[pos]--;
	hist[pointer[level1]] = pos;
	pointer[level1]++;	
	return;
}

void changePi(int level1, int *list, int *hist, int *pointer, int pos, int entry) {
	if (lambda==1) {
		list[pos]--;
		hist[pointer[level1]] = pos;
		pointer[level1]++;	

		Partial[level1-1]--;
	} else {
		list[pos] -= entry;
		hist[pointer[level1]] = pos;
		pointer[level1]++;	
		hist[pointer[level1]] = entry;
		pointer[level1]++;	
	
		Partial[level1-1] -= entry;
	}
	return;
}


void restore_entry_in_list(int level, int *list, int *hist, int *pointer, int sign) {
	int i;
	if (sign==-1) {
		for (i=pointer[level];i<pointer[level+1];i++) list[hist[i]]++;
	} else {
		for (i=pointer[level];i<pointer[level+1];i++) list[hist[i]]--;
	}
	pointer[level+1]=pointer[level];
	return;
}

void restore_Li(int level, int *list, int *hist, int *pointer) {
	int i;
	for (i=pointer[level];i<pointer[level+1];i++) list[hist[i]]++;
	pointer[level+1]=pointer[level];
	return;
}

void restorePi(int level, int *list, int *hist, int *pointer) {
	int i;
	if (lambda==1) {
		for (i=pointer[level];i<pointer[level+1];i++) {
			list[hist[i]]++;
			Partial[level]++;
		}
	} else {
		for (i=pointer[level];i<pointer[level+1];i+=2) {
			list[hist[i]] += hist[i+1];
			Partial[level] += hist[i+1];
		}
	}
	pointer[level+1]=pointer[level];
	return;
}

/******************************************************************
* The heart of the algorithm: the recursive function spread
*******************************************************************/
void spread(int level) {
	int i,j,p,l,k,li;
	int b;
	int numlines;
	int line;
	int point;
	
/******************************************************************
*  Some debugging output
*******************************************************************/
#if VERBOSE > 0
	printf("Level: %d\n",level);
/*	printf("Level: %d, jump_back: %d",level, jump_back);*/
#endif 

#if RANDOM	
	if (jump_back && level>0) return;
	if (jump_back && level==0) jump_back = 0;
#endif	

#if VERBOSE > 1
	for (i=0;i<numvar;i++) printf("%d ",Li[i]); printf("\n");
	for (i=0;i<numvar;i++) printf("%d ",Si[i]); printf("\n");
	for (i=0;i<numeqn;i++) printf("%d ",Pi[i]); printf("\n");
#endif 

/******************************************************************
*  Partial[level] and p contain the number of rows (n case lambda=1)
*  which are not handled yet.
*******************************************************************/
	if (level==0) {
		p = lambda*numeqn;
		Partial[level] = p;
		bestsofar = p;
	} else {
		Partial[level] = Partial[level-1];
		p = Partial[level];
	}

/******************************************************************
*  Print the partial solution vector
*******************************************************************/
#if VERBOSE > -1
if (loops%10000000==0) {
#if defined(FIX_FULL_ORBITS)
	printf("Loops: %ld, solutions: %ld, level: %d, bestsofar: %d, blockstodo: %ld\n",loops,nosolutions,level,bestsofar,blockstodo);
#else 	
	printf("Loops: %ld, solutions: %ld, level: %d, bestsofar: %d\n",loops,nosolutions,level,bestsofar);
#endif
	fflush(stdout);
}	
#endif

/******************************************************************
*  If p=0 then Pi is empty and Si is a solution vector.
*  Print the solution!
*******************************************************************/
	if (p<=0) {
#if defined(FIX_FULL_ORBITS)
		numberfulls = 0;
		for (i=0;i<numvar;i++) if (Si[i] && blocks[i]==order) numberfulls++; 
		if (numberfulls!=fullorbits) return;
#endif

	if (!SILENT) {
#if defined(SHORTOUTPUT)		
		for (i=0;i<numvar;i++) if (Si[i]==1) fprintf(solfile,"%d ",i); 
		fprintf(solfile,"\n");
#else
		for (i=0;i<numvar;i++) printf("%d",Si[i]); 
		printf(" L = %d\n",lambda); fflush(stdout);
		for (i=0;i<numvar;i++) fprintf(solfile,"%d",Si[i]); 
		fprintf(solfile,"\n");		
#endif		
		fflush(solfile);
}

/******************************************************************
*  We stop if we have reached the maximal number of solutions
*******************************************************************/
		nosolutions++;
		if ((stopafter>0)&&(stopafter<=nosolutions)) {
			printf("Loops: %ld\n",loops);
			printf("Stop after %ld solutions.\n",nosolutions);
			time_1 = os_ticks();
			user_time = time_1 - time_0;
			timestring[0] = 0;
			print_delta_time(user_time, timestring);
			printf("total enumeration time: %s\n", timestring);
			fflush(stdout);
			
			exit(1);
		}

		return;		 
	}

/******************************************************************
*  We stop if we have reached the maximal number of loops
*******************************************************************/
	loops++;
	if ((stoploops>0)&&(stoploops<loops)) {
		printf("Stop after %ld loops\n",loops);
		printf("Number of solutions: %ld\n",nosolutions);

		time_1 = os_ticks();
		user_time = time_1 - time_0;
		timestring[0] = 0;
		print_delta_time(user_time, timestring);
		printf("total enumeration time: %s\n", timestring);
		fflush(stdout);
		exit(1);
	}
/******************************************************************
*  Otherwise we have to go on.
*  There are still points in Pi to be covered, we have to go on.
*  Find a point p in Pi incident with the fewest - not yet 
*  excluded - lines Li(p) in Li. 
*******************************************************************/
	numlines = MAXINC;
	if (!switched) {
#if defined(FIX_FULL_ORBITS)
		for (i=0;i<numeqn;i++) {
			if (Pi[i]>0) { 
				b = 0; 
				for (j=Mrnum[i]-1;j>=0;j--) if ((Li[Mr[i][j]])&&(!Lifull[Mr[i][j]])) { 
					b++;
					if (b>=numlines) break;
				}
				if (b<numlines) {
					numlines = b;
					p = i;  	 
					if (numlines<=0) break;	 
				}	
			}
		}
#endif		
	} else {
		for (i=0;i<numeqn;i++) {
			if (Pi[i]>0) { 
				b = 0; 
				for (j=Mrnum[i]-1;j>=0;j--) if (Li[Mr[i][j]]) { 
					b++;
					if (b>=numlines) break;
				}
				if (b<numlines) {
					numlines = b;
					p = i;  	 
					if (numlines<=0) break;	 
				}	
			}
		}
	}
/*printf("p:%d numlines:%d %d\n",p,numlines,Pi[0]);*/

/******************************************************************
*  There is a point p in Pi incident with 0 lines Li(p) in Li. 
*  This means there is no solution for the partial solution
*  vector -> step back.
*******************************************************************/
	if (numlines==0) {

#if VERBOSE > 1
		printf("No line to choose -> return\n");
#endif
		return;
	}		

/******************************************************************
*  For the row (=point) p we store all incident lines, which are
*  not excluded yet, in Lp[level]
*******************************************************************/
	i = numlines-1;
	if (!switched) {
#if defined(FIX_FULL_ORBITS)
		for (l=Mrnum[p]-1;l>=0;l--) if ((Li[Mr[p][l]])&&(!Lifull[Mr[p][l]])) {
			Lp[level][i] = l;
			if (i==0) break;
			i--;
		}
#endif
	} else {
		for (l=Mrnum[p]-1;l>=0;l--) if (Li[Mr[p][l]]) {
			Lp[level][i] = l;
			if (i==0) break;
			i--;
		}
	}
	
#if VERBOSE > 1
	printf("Lines in L(%d): ",p);
	for (l=0;l<numlines;l++) printf("%d ",Lp[level][l]);
	printf("\n");
#endif	

/******************************************************************
*  For each l in Lp[level] we append the corresponding line to 
*  our partial solution and recurse.
*******************************************************************/
	for (l=0;l<numlines;l++) { 

#if RANDOM
again:			
		if (level < STARTBACKTRACK) l = rand() % numlines;
#endif		

		Pihistp[level+1] = Pihistp[level];
		Lihistp[level+1] = Lihistp[level];
		Sihistp[level+1] = Sihistp[level];

		line = Mr[p][Lp[level][l]];
#if VERBOSE > 0
		printf("Choose line: %d, Number of points: %d\n",line,Mcnum[line]);
		fflush(stdout);
#endif
/******************************************************************
* add line to the partial spread 
*******************************************************************/
		change_one_in_list(level+1,Si,Sihist,Sihistp,line,1); 
		
/******************************************************************
*  We treat the cases lambda=1 and lambda<>1 differently.
*  First, the case lambda=1
*******************************************************************/
		if (lambda==1) {
/******************************************************************
*  For each point in the line, we can mark the point in Pi
*  as done. 
*  Also, we all lines incident to this point cannot be contained in
*  a solution vector. This lines are excluded. We mark this in Li.
*******************************************************************/
			for (i=Mcnum[line]-1;i>=0;i--) {
				point = Mc[line][i];
				changePi(level+1,Pi,Pihist,Pihistp,point,1); 
				for (j=Mrnum[point]-1;j>=0;j--) {
					if (Li[Mr[point][j]])
/*					change_one_in_list(level+1,Li,Lihist,Lihistp,Mr[point][j],-1); */
						change_Li(level+1,Li,Lihist,Lihistp,Mr[point][j]); 
				}	
			}
		} else {
/******************************************************************
*  Second, the case lambda>1
*
*  We delete all previous chosen lines to avoid redundancy 
*  For each point in the line, we subtract the entry from Pi
*	We search for lines incident with the above point which are no
*  longer possible.
*******************************************************************/
			for (i=l;i>=0;i--) {
				k = Lp[level][i];
				if (Li[Mr[p][k]]==1) change_one_in_list(level+1,Li,Lihist,Lihistp,Mr[p][k],-1);
			}			
			for (i=Mcnum[line]-1;i>=0;i--) {
				point = Mc[line][i<<1];
				changePi(level+1,Pi,Pihist,Pihistp,point,Mc[line][(i<<1)+1]); 
			}
			for (i=Mcnum[line]-1;i>=0;i--) {
				point = Mc[line][i<<1];
				for (j=Mrnum[point]-1;j>=0;j--) {
					li = Mr[point][j];
					if (Li[li]) {
						for (k=Mcnum[li]-1;k>=0;k--) if (Mc[li][(k<<1)+1]>Pi[Mc[li][k<<1]]) {
							change_one_in_list(level+1,Li,Lihist,Lihistp,li,-1); 
							break;
						}
					}
				}	
			}
		}
		if (Partial[level]<bestsofar) bestsofar = Partial[level];
#if defined(FIX_FULL_ORBITS)
		blockstodo -= blocks[line];
		if (blockstodo==switchcapacity) {
			switched = 1;
			numberswitch++;
		}
#endif
/******************************************************************
*  Ok, we can step further on: Recursion
*******************************************************************/
 		spread(level+1);
/******************************************************************
*  The recursion step is done, we have to undo
*  our choices and go on to choose the next line of Lp
*******************************************************************/
		restorePi(level,Pi,Pihist,Pihistp);
/*		restore_entry_in_list(level,Li,Lihist,Lihistp,-1);*/
		restore_Li(level,Li,Lihist,Lihistp);
/*
		if (!switched)
			restore_Li(level,Lifull,Lifullhist,Lifullhistp);
*/			
		restore_entry_in_list(level,Si,Sihist,Sihistp,1);
#if defined(FIX_FULL_ORBITS)
		if (blockstodo==switchcapacity) {
			switched = 0;
		}		
		blockstodo += blocks[line];
#endif

#if RANDOM
		if (level <= STARTBACKTRACK && level > 0) {
			jump_back = 1;
			return;
		}
		if (level == 0) { jump_back = 0; goto again; }
#endif		
	}
/******************************************************************
*  If all lines are done we step back. 
*******************************************************************/
	return;
}


/******************************************************************
*  The main program
*******************************************************************/
int main(int argc, char *argv[]) {
	int i,j,k,s,b,i_old;
	int b1;
	int j0;
	int level;
	char zeile[zlength];
	char *filename;
	char select[100];
	char *zp;
	
/******************************************************************
*  The command line parameters are read
*******************************************************************/
	if ((argc<3)&&(argc>5)) {
		printf("Wrong parameters! Correct are: lambda [silent] filename [SEL]\n");
		exit(0) ;
	}
 	lambda = atoi(argv[1]);

  selection = SILENT = 0;
  for (i=1;i<argc;i++) {
		if (strcmp(argv[i],"silent")==0) SILENT = 1;
  }
	if (argc == 4+SILENT)
		selection = atoi(argv[argc-1]);

	filename = argv[2+SILENT];
	nosolutions = 0;
	loops = 0;

/******************************************************************
* Look for maximal number of solutions or loops
*******************************************************************/
	stopafter = -1;
	stoploops = -1;
#if defined(FIX_FULL_ORBITS)
	fullorbits = 0;
#endif
	infile = fopen(filename,"r");
	do {
		fgets(zeile,zlength,infile);
      if (strstr(zeile,"% stopafter")!=NULL) {
			sscanf(zeile,"%% stopafter %ld",&stopafter);
		}
      if (strstr(zeile,"% stoploops")!=NULL) {
			sscanf(zeile,"%% stoploops %ld",&stoploops);
		}
#if defined(FIX_FULL_ORBITS)
      if (strstr(zeile,"% fullorbits")!=NULL) {
			sscanf(zeile,"%% fullorbits %d",&fullorbits);
		}
      if (strstr(zeile,"% order:")!=NULL) {
			sscanf(zeile,"%% order: %ld",&order);
		}
#endif
	} while (zeile[0]=='%');
	sscanf(zeile,"%d%d",&numeqn,&numvar);

/******************************************************************
* Look for number of equations and variables
*******************************************************************/
	printf("numvar: %d, numeqn: %d\n",numvar,numeqn);
	maxdepth = 2 + (numeqn*lambda<numvar)?numeqn*lambda:numvar;

#if defined(FIX_FULL_ORBITS)
	printf("fullorbits: %d\n",fullorbits);
#endif
		
/******************************************************************
* Allocate memory 
*******************************************************************/
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
#if defined(FIX_FULL_ORBITS)
	Lifull = (int*)calloc(numvar,sizeof(int));
	blocks = (long*)calloc(numvar,sizeof(long));
#endif

	Partial = (int*)calloc(maxdepth,sizeof(int));
	
	Mcnum = calloc(numvar,sizeof(int));
	Mrnum = calloc(numeqn,sizeof(int));

	for (i=0;i<numvar;i++) Si[i]=1;
	for (i=0;i<numvar;i++) Mcnum[i] = 0;
	for (i=0;i<numeqn;i++) Mrnum[i] = 0;

/******************************************************************
* Read the matrix and write the nonzero entries to 
* spreadtmp.txt. All columns with entries > lambda are deleted.
* Also the number of entries for each line and each column
* are counted in order to save space.
*******************************************************************/
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

#if defined(FIX_FULL_ORBITS)
	for (j=0;j<numvar;j++) {
		fscanf(infile,"%ld",&(blocks[j]));
		if (blocks[j]==order) Lifull[j] = 1;
	}
	fscanf(infile,"%ld",&numberofblocks);
#endif
/******************************************************************
* The selection vector is handled
*******************************************************************/
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

/******************************************************************
* Now, the nonzero entries have been written into the file
* spreadtmp.txt. We allocate memory for Mr and Mc.
*******************************************************************/
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

	Mr = calloc(numeqn,sizeof(int*));
	for (i=0;i<numeqn;i++) {
		Mr[i] = calloc(Mrnum[i],sizeof(int));
	}

	Mc = calloc(numvar,sizeof(int*));
	if (lambda==1) {
		for (i=0;i<numvar;i++) {
			Mc[i] = calloc(Mcnum[i],sizeof(int));
		}
	} else {
		for (i=0;i<numvar;i++) {
			Mc[i] = calloc(2*Mcnum[i],sizeof(int));
		}
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

/******************************************************************
*  The entries from are read from spreadtmp.txt
*******************************************************************/
	for (i=0;i<numvar;i++) Li[i] = 0;
	j0 = 0;
	printf("deleting unnecessary lines... \n");
	fflush(stdout);	
	for (k=0;k<b;k++) {
		i_old = i;
		fscanf(infile,"%d%d%d",&i,&j,&s);
		if (i!=i_old) j0 = 0;
		if (Si[j]==1) {		
			Mr[i][j0] = j;
			j0++;
			
			if (lambda==1) {			
				Mc[j][Li[j]] = i;
			} else {
				Mc[j][(Li[j]<<1)+1] = s;
				Mc[j][Li[j]<<1] = i;
			}
			Li[j]++;
		} else {
			Mcnum[j] = 0;
			Mrnum[i]--;
		}
	}
	fclose(infile);

	printf("Input finished\n");
	fflush(stdout);
	
/******************************************************************
* Now the matrix has been read again and it is available twice:
* via column number and via row number in Mc resp. Mr.
*
* Initially, the arrays for points and lines are cleared. 
* The partial spread is empty. 
*******************************************************************/
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
#if defined(FIX_FULL_ORBITS)
	for (i=0;i<numvar;i++) if (!Li[i]) Lifull[i]=0;
/*	for (i=0;i<numvar;i++) if (Lifull[i]) Li[i]=0;*/
	switchcapacity = fullorbits*order;
	blockstodo = numberofblocks;
	switched = 0;
	printf("order: %ld numberofblocks: %ld switchcapacity: %ld\n",order,numberofblocks,switchcapacity);
	numberswitch=0;
#endif

#if VERBOSE > 0	
	printf("Start of recursion\n"); fflush(stdout);
#endif	
	if (SILENT) {
		printf("No output of solutions !!!\n"); fflush(stdout);
	} else {
		printf("Output of solutions into file solutions\n"); fflush(stdout);
		solfile = fopen("solutions","w");
	}

/******************************************************************
*  Ok, we can start the recursive loop
*******************************************************************/
	time_0 = os_ticks();
	spread(level);
	time_1 = os_ticks();

/******************************************************************
*  Now we are done. We print some statistics and exit
*******************************************************************/	
	if (!SILENT) fclose(solfile);
	printf("Loops: %ld\n",loops);
	printf("\nNumber of solutions: %ld\n",nosolutions);
	
	user_time = time_1 - time_0;
	timestring[0] = 0;
	print_delta_time(user_time, timestring);
	printf("total enumeration time: %s\n", timestring);
	fflush(stdout);

	return 1;
}
