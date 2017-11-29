// args.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include "packcover.H"
#include "args.H"
#include "utils.H"
#include "warning.H"

OptionLine::OptionLine(int argc, char **argv)
{
  if(!(argc & 1)) errorExit("odd number of arguments");
  for(int i = 1; i < argc; i += 2) {
    args[argv[i]].push_back(string(argv[i+1]));
    used[string(argv[i])] = false;
  }
}

OptionLine::~OptionLine()
{
  for(map<string,bool>::const_iterator I = used.begin(); I != used.end(); ++I)
    if(!I->second)
      warning->print(string("Unused option: ") + I->first);
}
