#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "dbms.h"

FILE *log;
FILE *LogFile; /* Aenderung wolfb: zusaetzliche Zeile! */

int main( int argc, char *argv[] ) {
   float  t0,t1;
   char  *l;
   int    error;

   if ( argc < 3 ) {
      printf("DBREMOVE - Deleting records\n");
      printf("Call: dbremove <record file> <table>\n");
   } else {
      printf("MOLGEN DataBase Version 1.0   Juergen Mueller 1993\n");
      t0 = clock();
      error = dbms_remove_f( argv[2], argv[1], 1 );
      if ( error ) { printf("Operation failed, return code %d\n", error ); }
      t1 = clock();
      printf("Time required: %.0f seconds\n", (t1-t0) / CLOCKS_PER_SEC /* CLK_TCK */ );
   }
}
