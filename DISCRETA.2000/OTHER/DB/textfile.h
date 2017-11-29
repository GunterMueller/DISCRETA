#define TFBUFSIZE 1024

#define NO_LEADING_SPACES 1
#define NO_ENDING_SPACES  2
#define NO_LEADEND_SPACES 3

class TextFile {

   private:

      FILE     *file;
      unsigned  bufsize;
      char     *buffer1;
      char     *buffer2;

   public:

      TextFile();
     ~TextFile();

      char *init( char *name );
      char *init( char *name, unsigned size );
      char *readline();
      char *readline( char *buf );
      char *writeline();
      char *skiplines( unsigned cnt );
      char *getfield( unsigned pos, unsigned len, unsigned mode );
      void  drop();
};
