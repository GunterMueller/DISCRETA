/* constructing error-correcting codes with a nontrivial 
   permutation group;
   transforms problem into a clique instance,   
   by Patric Ostergard, 30.3.1999 */

/* gcc group2graph.c -o group2graph -O2 */

#define INT_SIZE (8*sizeof(unsigned))
#define MAX_N 30
#define SPACE_SIZE 10000000  /* > 2^23 => works for n <= 28 */
#define MAX_PERM 10000
#define MAX_TAB 50000
#define marked(a) (bit[a/INT_SIZE]&(mask[a%INT_SIZE]))

#include <stdio.h>

int n,w,d,groupsize; /* parameters, q=2 !!! */
int perm[MAX_PERM][MAX_N];
int table[MAX_TAB][MAX_N],orbsize[MAX_TAB];
int exp[MAX_N];
int binCoef[MAX_N][MAX_N+1];
unsigned mask[INT_SIZE],bit[SPACE_SIZE];

mark(a)
int a;
{
    bit[a/INT_SIZE]|=mask[a%INT_SIZE];
}

/* rank and unrank algorithms */

firstset(rep,all)
int rep[MAX_N],all;
{
  int i;

  if(all) {
    for(i=0;i<n;i++)
      rep[i] = 0;
  }
  else {
    for(i=0;i<w;i++)
      rep[i] = 1;
    for(i=w;i<n;i++)
      rep[i] = 0;
  }
}

nextset(rep,all)
int rep[MAX_N],all;
{   
  int i,j,k,pos,rank;

  do {
    if(all) {
      pos = n-1;
      do {
        rep[pos] = (rep[pos]+1)%2;
        pos--;
      } while((rep[pos+1] == 0) && (pos!=-1));
      if((pos==-1)&&(rep[pos+1]==0)) return 0;
    }
    else {
      i = 0;
      while(rep[i]==0)
        i++;  
      if(i>=n-w) return 0;
      else {
        j = i+1;
        while(rep[j]==1)
          j++;
        rep[j] = 1;
        for(k=0;k<j-i-1;k++)
          rep[k] = 1;
        for(k=j-i-1;k<j;k++)
          rep[k] = 0;
      }
    }
    rank = 0;
    for(i=0;i<n;i++)
      rank += exp[i]*rep[i];
  } while(marked(rank));
  return 1;
}

int minlex(wd)
int wd[MAX_N];
{
  int i,j;
  int rank,rank2,orb_size;
  int orbel[MAX_PERM],dist;
  int reject;

  reject = 0;
  orb_size = 1;
  rank = 0;
  for(i=0;i<n;i++)
    rank += exp[i]*wd[i];
  orbel[0] = rank;
  mark(rank);
  for(i=1;i<groupsize;i++) {
    rank2 = 0;
    for(j=0;j<n;j++)
      rank2 += exp[j]*wd[perm[i][j]];
    dist = 0;
    for(j=0;j<n;j++)
      if(wd[perm[i][j]]!=wd[j]) dist++;
    if((dist<d)&&(dist>0)) reject = 1; /* min. dist. too small */
    else {    
      for(j=0;j<orb_size;j++)
        if(orbel[j]==rank2) goto end;
      orbel[orb_size++] = rank2;
    }
    mark(rank2);
end: ;
  }
  if(reject) return 0;
  else return orb_size;
}

int dist(wd1,wd2)
int wd1[MAX_N],wd2[MAX_N];
{
  int i,j;
  int dst,mindist;

  mindist = 1000;
  for(i=0;i<groupsize;i++) {
    dst = 0;
    for(j=0;j<n;j++)
      if(wd1[perm[i][j]]!=wd2[j]) dst++;
    if(dst<mindist) mindist = dst;
  }
  return mindist;
}

int mult(a,b,c,d)
int a,b[MAX_N],c[MAX_N],d[MAX_N];
{
  int i;
  
  for(i=0;i<a;i++)
    d[i] = b[c[i]];
}

main(argc,argv)
int argc;
char *argv[];
{
  register int i,j,k,k2;
  int count,count2,temp,pos,quasi,temp1,temp2;
  int rep[MAX_N],temptab[MAX_TAB];
  int new[MAX_PERM][MAX_N],g[MAX_PERM][MAX_N]; 
  int last[MAX_PERM][MAX_N],gen[MAX_PERM][MAX_N];
  int f[MAX_N];
  int new_ant,g_ant,last_ant,gen_ant,all,glength;

  if(argc!=4) 
    err_message();
  if((n = atoi(argv[1]))==0)
    err_message();
  if((d = atoi(argv[2]))==0)
    err_message();
  if((w = atoi(argv[3]))==0)
    err_message();
  if(w==1) all = 1;
  else all = 0;

  /* initialize mask */
  mask[0] = 1;
  for(i=1;i<INT_SIZE;i++)
    mask[i] = mask[i-1]<<1;

  /* initialize tables */
  for(i=0;i<SPACE_SIZE;i++)
    bit[i] = 0;    

  for(i=0;i<MAX_N;i++) {
    binCoef[i][0] = 1;
    binCoef[i][i] = 1;
    binCoef[i][i+1] = 0;
    for(j=1;j<i;j++)
      binCoef[i][j] = binCoef[i-1][j-1]+binCoef[i-1][j];
  }
  exp[0] = 1;
  for(i=1;i<n;i++)
    exp[i] = exp[i-1]*2;

  /* read in generators */
  if(!scanf("%d %d\n",&gen_ant,&glength))
    err_message(); 
  for(i=0;i<gen_ant;i++)
    for(j=0;j<glength;j++) {
      if(!scanf("%d",&temp)) 
        err_message();
      gen[i][j] = temp-1; /* DISCRETA: 1->  Here: 0-> */
    }

  /* this is Algorithm 6.5 in Kreher & Stinson */
  g_ant = 0;
  new_ant = 1;
  for(i=0;i<n;i++)
    new[0][i] = i;
  while (new_ant > 0) {
    for(i=0;i<new_ant;i++)
      for(j=0;j<n;j++)
        g[g_ant+i][j] = new[i][j];
    g_ant += new_ant;
    for(i=0;i<new_ant;i++)
      for(j=0;j<n;j++)
        last[i][j] = new[i][j];
    last_ant = new_ant;
    new_ant = 0;
    for(i=0;i<gen_ant;i++)
      for(j=0;j<last_ant;j++) {
        mult(n,&gen[i][0],&last[j][0],f);
        for(k=0;k<g_ant;k++) {
          for(k2=0;k2<n;k2++)
            if(f[k2]!=g[k][k2]) goto con;
          goto exists;
      con: ;
	}
	/* also check whether it is in new */
        for(k=0;k<new_ant;k++) {
          for(k2=0;k2<n;k2++)
            if(f[k2]!=new[k][k2]) goto con2;
          goto exists;
      con2: ;
	}
        for(k=0;k<n;k++)
          new[new_ant][k] = f[k];
        new_ant++;
  exists: ;       
      }
  }
  for(i=0;i<g_ant;i++) {
    for(j=0;j<n;j++) {
      perm[i][j] = g[i][j];
/*    printf("%d ",g[i][j]); */
    }
/*  printf("\n"); */
  }
  groupsize = g_ant;
  /* printf("%d\n",groupsize);
  exit(); */

  /* Algorithm 6.5 ends */
  /* initialization ready */
  /* start main routine */

  firstset(rep,all);
  count = 0;
  j = 0;
  do {
    temp = minlex(rep); /* also checks internal distance */
    if(temp==0) goto next;
    for(i=0;i<n;i++)
      table[count][i] = rep[i];
    orbsize[count++] = temp;
next:;
  } while(nextset(rep,all));
  printf("%d 1000\n",count);
  for(i=0;i<count;i++) {
    count2 = 0;
    for(j=0;j<count;j++) {
      if(i==j) continue;
      if(dist(&table[i][0],&table[j][0])>=d) temptab[count2++] = j;
    }
    printf("%d %d ",orbsize[i],count2);
    for(j=0;j<count2;j++)
      printf("%d ",temptab[j]);
    printf("\n");
  }
}

err_message()
{
printf("Usage: group2graph n d w\n"); 
exit(1);
}


