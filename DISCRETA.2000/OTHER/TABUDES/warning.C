// warning.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include "packcover.H"
#include "warning.H"

void Warning::print(const string& msg)
{
  int now = ++soFar[msg];
  if(now < printCount)
    out << "Warning: " << msg << '\n';
  else if(now == printCount)
    out << "Warning: " << msg << "\n(suppressed until exit of program)\n";
}

void Warning::flush()
{
  if(soFar.empty()) return;
  out << "Summary of warnings:\n";
  for(map<string,int>::const_iterator I = soFar.begin();
      I != soFar.end(); ++I)
    out << I->second << " warning" << (I->second > 1 ? "s" : "") << ": "
	<< I->first << '\n';
}
