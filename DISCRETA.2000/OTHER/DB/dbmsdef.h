/*--------------------------------------------------------------------------*/
/*                                                                          */
/* DBMSDEF.H - Definition wichtiger Datenbankparameter                      */
/*                                                                          */
/*--------------------------------------------------------------------------*/

#include "dbmserr.h"

extern FILE *LogFile;
void Log(char *s);
void *Memmove(void *s1, const void *s2, int len);
void Memory_zero(void *p, int l);

/* Directory separator */
/* DOS or OS/2 #define DIRSEP                   '\\' */
/* UNIX        #define DIRSEP                   '/' */

#if 0
/* AB 950118 */
#define DIRSEP                   '\\'
#endif
#define DIRSEP                   '/'

/* Extensions of DBMS Files */
#define TBLEXT                   ".tbl"
#define TDFEXT                   ".tdf"
#define DEFEXT                   ".def"

#define FNSIZE                   256
#define TBLINF                   80
#define BASCNT                   10
#define INDEXC                   32
#define BASE_PAGE_SIZE           1024
#define FILE_BEGIN               SEEK_SET
#define DP  1
#define OP  2
#define IP  3

/*--------------------------------------------------------------------------*/
/*                                                                          */
/* Defintion of table operations                                            */
/*                                                                          */
/*--------------------------------------------------------------------------*/

#define FIO         1
#define INSERT      2

/*--------------------------------------------------------------------------*/
/*                                                                          */
/* Definition of index specials                                             */
/*                                                                          */
/*--------------------------------------------------------------------------*/

#define AUTOINCREMENT 1
#define AUTODECREMENT 2


/*--------------------------------------------------------------------------*/
/*                                                                          */
/* Defintion of basic types                                                 */
/*                                                                          */
/*--------------------------------------------------------------------------*/

#define CHAR   char
#define SCHAR  signed char
#define UCHAR  unsigned char
#define SINT   signed int
#define UINT   unsigned int
#if 0
#define SLONG  signed long
#define ULONG  unsigned long
#endif

/* AB 950118 */
#define SLONG  signed int
#define ULONG  unsigned int

#define FLOAT  float
#define DOUBLE double
#define PSIZE  unsigned int


#define DB_ALL_OK                0

#define DB_EXIST_ERROR           1

#define DB_MEMORY_ERROR          8
#define DB_ACCESS_DENIED         9
#define DB_DATA_TOO_LARGE       10
#define DB_STACK_ERROR          11
#define DB_UNDEFINED_CLASS      12
#define DB_INPUT_ERROR          13
#define DB_OUTPUT_ERROR         14
#define DB_PAGE_SIZE_ERROR      17
#define DB_DATA_NOT_FOUND       18
#define DB_WRONG_FILENAME       19
#define DB_SPLIT_PAGE           100
#define DB_SPLIT_CRASH          111

#define DB_INDEX_PAGE 0
#define DB_DATA_PAGE 1

#define MIN_PAGE_SIZE 0x0100
#define MAX_PAGE_SIZE 0xFFFF


#define LEFT  1
#define RIGHT 0

#if 0
#define DOS
#endif
/* AB 950118 */

#define _GREATER 0
#define _SMALLER 1


#define _CREATE_
#define _INSERT_
#define _CHANGE_
#define _DELETE_
#define _SEARCH_
#define _VIEW_

