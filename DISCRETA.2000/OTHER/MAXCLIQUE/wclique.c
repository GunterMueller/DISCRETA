/* wclique.c exact algorithm for finding one maximum-weight 
   clique in an arbitrary graph,
   21.1.1999, Patric R. J. Ostergard, 
   patric.ostergard@hut.fi */

/* compile: gcc wclique.c -o wclique -O2 */

/* usage: wclique infile */

/* infile format: see http://www.tcs.hut.fi/~pat/wclique.html */

#include <stdio.h>
#include <sys/times.h>
#include <sys/types.h>

#define INT_SIZE (8*sizeof(int))
#define TRUE 1
#define FALSE 0
#define MAX_VERTEX 7000 /* maximum number of vertices */
#define is_edge(a,b) (bit[a][b/INT_SIZE]&(mask[b%INT_SIZE]))

int Vnbr,Enbr;          /* number of vertices/edges */
int clique[MAX_VERTEX]; /* table for pruning */
int bit[MAX_VERTEX][MAX_VERTEX/INT_SIZE+1];
int wt[MAX_VERTEX];

int pos[MAX_VERTEX];    /* reordering function */
int set[MAX_VERTEX];    /* current clique */
int rec[MAX_VERTEX];	/* best clique so far */
int record;		/* weight of best clique */
int rec_level;          /* # of vertices in best clique */

unsigned mask[INT_SIZE];
void graph();           /* reads graph */
void write_result();    /*  */

struct tms bf;
int timer1;
double timer11;

main (argc,argv)
int argc;
char *argv[];
{
  int i,j,k,p;
  int max_wt,max_nwt,wth;
  int new[MAX_VERTEX],used[MAX_VERTEX];
  int nwt[MAX_VERTEX];
  FILE *infile;

  /* read input */
  if(argc < 2) {   	
    printf("Usage: wclique infile\n");
    exit(1);
  }
  if((infile=fopen(argv[1],"r"))==NULL)
    fileerror();    

  /* initialize mask */
  mask[0] = 1;
  for(i=1;i<INT_SIZE;i++)
    mask[i] = mask[i-1]<<1;

  /* read graph */
  graph(infile);
 
  /* "start clock" */
  times(&bf);
  timer1 = bf.tms_utime;

  /* order vertices */
  for(i=0;i<Vnbr;i++) {
    nwt[i] = 0;
    for(j=0;j<Vnbr;j++)
      if (is_edge(i,j)) nwt[i] += wt[j];
  }
  for(i=0;i<Vnbr;i++)
    used[i] = FALSE;
  for(i=Vnbr-1;i>=0;i--) {
    max_wt = -1; max_nwt = -1;
    for(j=0;j<Vnbr;j++)
      if((!used[j])&&((wt[j]>max_wt)||(wt[j]==max_wt&&nwt[j]>max_nwt))) {
        max_wt = wt[j];
        max_nwt = nwt[j];
        p = j;
      }
    pos[i] = p;
    used[p] = TRUE;
    for(j=0;j<Vnbr;j++)
      if((!used[j])&&(j!=p)&&(is_edge(p,j)))
        nwt[j] -= wt[p];
  }

  /* main routine */
  record = 0;
  wth = 0;
  for(i=0;i<Vnbr;i++) {
     wth += wt[pos[i]];
     sub(i,pos,0,0,wth);
     clique[pos[i]] = record;
     times(&bf);
     timer11 = (bf.tms_utime - timer1)/100.0;
     printf("level = %3d(%d) best = %2d time = %8.2f\n",i+1,Vnbr,record,timer11);
	  fflush(stdout);
  }
  printf("Record: ");
  max_wt = 0;
  for(i=0;i<rec_level;i++) {
    printf ("%d ",rec[i]);
	 max_wt += wt[rec[i]];
  }
  printf ("\n");
  printf ("Blocks: %d\n",max_wt);
  write_result();
}

int sub(ct,table,level,weight,l_weight)
int ct,level,weight,l_weight;
int *table;
{
  register int i,j,k;
  int best;
  int curr_weight,left_weight;
  int newtable[MAX_VERTEX];
  int *p1,*p2;

  if(ct<=0) { /* 0 or 1 elements left; include these */
    if(ct==0) { 
      set[level++] = table[0];
      weight += l_weight;
    }
    if(weight>record) {
      record = weight;
      rec_level = level;
      for (i=0;i<level;i++) rec[i] = set[i];
    }
    return 0;
  }
  for(i=ct;i>=0;i--) {
    if((level==0)&&(i<ct)) return 0;
    k = table[i];
    if((level>0)&&(clique[k]<=(record-weight))) return 0;  /* prune */
    set[level] = k;
    curr_weight = weight+wt[k];
    l_weight -= wt[k];
    if(l_weight<=(record-curr_weight)) return 0; /* prune */
    p1 = newtable;
    p2 = table;
    left_weight = 0;   
    while (p2<table+i) { 
      j = *p2++;
      if(is_edge(j,k)) {
		  *p1++ = j;
        left_weight += wt[j];
      }
    }
    if(left_weight<=(record-curr_weight)) continue;
    sub(p1-newtable-1,newtable,level+1,curr_weight,left_weight);
  }
  return 0;
}

void graph(fp) 
FILE *fp;
{
  register int i,j,k;
  int weight,degree,entry;
 
  if(!fscanf(fp,"%d %d\n",&Vnbr,&Enbr))
    fileerror(); 
  for(i=0;i<Vnbr;i++)     /* empty graph table */
    for(j=0;j<Vnbr/INT_SIZE+1;j++)
      bit[i][j] = 0;
  for(i=0;i<Vnbr;i++) {
    if(!fscanf(fp,"%d %d",&weight,&degree))
      fileerror(); 
    wt[i] = weight;
    for(j=0;j<degree;j++) {
      if(!fscanf(fp,"%d",&entry))
        fileerror();
      bit[i][entry/INT_SIZE] |= mask[entry%INT_SIZE]; /* record edge */
    }
  }
  fclose(fp);
}

void write_result() {
	FILE *outfile;
	int i;

	outfile = fopen("outfile","w");
	fprintf (outfile,"%d\n",rec_level);
	for(i=0;i<rec_level;i++) {
		fprintf (outfile,"%d ",rec[i]);
	}
	fclose(outfile);
	return;
}

int fileerror()
{
  printf("Error in graph file\n");
  exit();
}
