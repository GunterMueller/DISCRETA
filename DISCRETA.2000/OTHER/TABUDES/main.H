// main.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_main_h
#define PACKCOVER_main_h

#include <stack.h>
#include <string>
#include "utils.H"

class NeighParseStackElem {
private:
  enum neighStackType { nptr, str };
  const neighStackType typeId;
  Neighborhood* neighPtr;
  string tokenStr;
public:
  NeighParseStackElem(Neighborhood *ptr) : typeId(nptr), neighPtr(ptr) {}
  NeighParseStackElem(string s) : typeId(str), tokenStr(s) {}
  void parseError() { errorExit("neighborhoods: parse error"); }
  string getStr() { 
    if(typeId != str) parseError();
    return tokenStr;
  }
  Neighborhood* getNptr() {
    if(typeId != nptr) parseError();
    return neighPtr;
  }
};

class AttrParseStackElem {
private:
  enum attrStackType { aptr, str };
  const attrStackType typeId;
  Attribute* attrPtr;
  string tokenStr;
public:
  AttrParseStackElem(Attribute *attr) : typeId(aptr), attrPtr(attr) {}
  AttrParseStackElem(string s) : typeId(str), tokenStr(s) {}
  void parseError() { errorExit("attributes: parse error"); }
  string getStr() { 
    if(typeId != str) parseError();
    return tokenStr;
  }
  Attribute* getAptr() {
    if(typeId != aptr) parseError();
    return attrPtr;
  }
};

#endif
