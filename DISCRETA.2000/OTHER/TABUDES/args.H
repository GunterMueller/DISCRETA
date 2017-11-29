// args.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_args_h
#define PACKCOVER_args_h

#include <string>
#include <strstream.h>
#include <map.h>
#include <list.h>

class OptionLine {
private:
  map<string,list<string> > args;
  map<string,bool> used;
public:
  OptionLine(int argc, char **argv);
  ~OptionLine();

  // getNext returns the next value of that name, or returns false if
  // such value is not found. getNextMust returns the next value of
  // that name, and aborts if no such value is defined

  template<class T>
  bool getNext(const char *name, T& value) {
    used[name] = true;
    if(args[string(name)].empty())
      return false;
    string val(args[string(name)].front());
    istrstream istr(val.c_str(), val.length());
    if(!(istr >> value))
      errorExit(string("unexpected format for argument: ") +
		string(name));
    args[string(name)].pop_front();
    return true;
  }

  template<>
  bool getNext(const char *name, string& value) {
    used[name] = true;
    if(args[string(name)].empty())
      return false;
    value = args[string(name)].front();
    args[string(name)].pop_front();
    return true;
  }

  template<class T>
  T getNextMust(const char *name, const T dummy) {
    used[name] = true;
    T tmp;
    if(!getNext(name, tmp)) 
      errorExit(string("Required option <") + string(name) + 
		string("> not found"));
    return tmp;
  }

  // get returns the first value of that name, or returns default, if
  // value is not defined

  template<class T>
  T get(const char *name, const T deflt) {
    used[name] = true;
    if(args[string(name)].empty())
      return deflt;
    string val(args[string(name)].front());
    istrstream istr(val.c_str(), val.length());
    T tmp;
    if(!(istr >> tmp))
      errorExit(name);
    return tmp;
  }

  template<>
  string get(const char *name, const string deflt) { 
    used[name] = true;
    if(args[string(name)].empty())
      return deflt;
    else
      return args[string(name)].front();
  }

  template<class T>
  T getMust(const char *name, const T dummy) {
    used[name] = true;
    if(args[string(name)].empty())
      errorExit(string("Required option <") + string(name) + 
		string("> not found"));
    return get(name, dummy);
  }
};

#endif
