/* kmgen.c */

/*  Brute force Kramer-Mesner matrix generation.
**
**  (c) 1997,1999 Kari J. Nurmela
**
** usage: kmgen permutation-file order k t Kramer-Mesner-file
*/

#include <stdio.h>

int main(int argc, char *argv[])
{
  FILE *pfile;
  int permcount, pointcount;
  int i;
  char buf[1000], *ptr;

  if(argc != 6)
    abort();
  pfile = fopen(argv[1], "r");
  fscanf(pfile, "%d %d\n", &permcount, &pointcount);
  printf("Pack %d %s %s ", pointcount, argv[3], argv[4]);
  for(i = 0; i < permcount; ++i) {
    ptr = buf;
    while(1) {
      fscanf(pfile, "%c", ptr);
      if(*ptr == '\n') {
	*ptr = '\0';
	break;
      }
      else
	++ptr;
    }
    printf("\"pres %s\" ", buf);
  }
  printf("-o -pt -ns -kmfile %s %s\n", argv[5], argv[2]);
  return 0;
}
