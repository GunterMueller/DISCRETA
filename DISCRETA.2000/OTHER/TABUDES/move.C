// move.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include "packcover.H"
#include "move.H"
#include "utils.H"
#include "solution.H"

#include <iostream.h>

ostream& operator<<(ostream& out, Move& move)
  // mainly for general neighborhood debug purposes
{
  out << "Move: add ";
  copy(move.getAdded().begin(), move.getAdded().end(), 
       ostream_iterator<patchType>(out, " "));
  out << "remove ";
  copy(move.getRemoveIndexes().begin(), move.getRemoveIndexes().end(), 
       ostream_iterator<int>(out, " "));
  return out;
}

void Move::print(Solution& solution, ostream& out)
{
  out << "Move: add patches ";
  copy(getAdded().begin(), getAdded().end(), 
       ostream_iterator<patchType>(out, " "));
  out << "remove patches ";
  for(Move::c_indexIterator I = getRemoveIndexes().begin();
      I != getRemoveIndexes().end(); ++I)
    out << solution.getPatch(*I) << ' ';
}

void Move::append(const Move& move)
{
  insert_iterator<list<patchType> > patchInsert(added, added.end());
  copy(move.getAdded().begin(), move.getAdded().end(), patchInsert);
  addCnt += move.addCnt;
  insert_iterator<list<int> > indexInsert(removed, removed.end());
  copy(move.getRemoveIndexes().begin(), move.getRemoveIndexes().end(),
       indexInsert);
  remCnt += move.remCnt;
}

void Move::checkConsistency()
{
  // check for double removals
  if(removed.size() < 2) return;
  vector<int> tmp(removed.size());
  copy(removed.begin(), removed.end(), tmp.begin());
  sort(tmp.begin(), tmp.end());
  if(adjacent_find(tmp.begin(), tmp.end()) != tmp.end())
    errorExit("move with duplicate removal indices");
}

void Move::printRemovePatches(ostream& out, Solution& solution)
{
  for(Move::c_indexIterator I = getRemoveIndexes().begin();
      I != getRemoveIndexes().end();++I)
    out << solution.getPatch(*I) << " ";
}
