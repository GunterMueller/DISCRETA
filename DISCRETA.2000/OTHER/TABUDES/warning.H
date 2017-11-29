// warning.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_warning_h
#define PACKCOVER_warning_h

#include <map.h>
#include <string>

class Warning {
private:
  ostream& out;
  const int printCount;
  map<string,int> soFar;
public:
  Warning(ostream& os, int pc = 4) : out(os), printCount(pc) {}
  void print(const string& msg);
  void flush();
};

#endif
