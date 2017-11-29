// attr.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include "packcover.H"
#include "attr.H"
#include "move.H"
#include "solution.H"
#include "neighborhood.H"

int AttrIdProvider::attrIds = 0;

list<vector<int> > RemoveAttribute::makeAttr
(Move& m, Solution&, volumeType, penaltyType) const
{
  list<vector<int> > retval;
  for(Move::c_patchIterator I=m.getAdded().begin(); 
      I != m.getAdded().end(); ++I)
    retval.push_back(makeVector(int(*I), getAttrId()));
  return retval;
}

list<vector<int> > RemoveAttribute::makeTestAttr
(Move& m, Solution& solution, volumeType, penaltyType) const
{
  list<vector<int> > retval;
  for(Move::c_indexIterator I=m.getRemoveIndexes().begin();
      I != m.getRemoveIndexes().end(); ++I)
    retval.push_back(makeVector(int(solution.getPatch(*I)), getAttrId()));
  return retval;
}

list<vector<int> > AddAttribute::makeAttr
(Move& m, Solution& solution, volumeType, penaltyType) const
{
  list<vector<int> > retval;
  for(Move::c_indexIterator I=m.getRemoveIndexes().begin();
      I != m.getRemoveIndexes().end(); ++I)
    retval.push_back(makeVector(int(solution.getPatch(*I)), getAttrId()));
  return retval;
}

list<vector<int> > AddAttribute::makeTestAttr
(Move& m, Solution&, volumeType, penaltyType) const
{
  list<vector<int> > retval;
  for(Move::c_patchIterator I=m.getAdded().begin();
      I != m.getAdded().end(); ++I)
    retval.push_back(makeVector(int(*I), getAttrId()));
  return retval;
}

list<vector<int> > AddAllAttribute::makeAttr
(Move& m, Solution& solution, volumeType, penaltyType) const
{
  list<vector<int> > retval;
  vector<int> removed(m.removeCount() + 1);
  vector<int>::iterator vI = removed.begin();
  for(Move::c_indexIterator I=m.getRemoveIndexes().begin();
      I != m.getRemoveIndexes().end(); ++I, ++vI)
    *vI = solution.getPatch(*I);
  *vI = getAttrId();
  retval.push_back(removed);
  return retval;
}

list<vector<int> > AddAllAttribute::makeTestAttr
(Move& m, Solution&, volumeType, penaltyType) const
{
  list<vector<int> > retval;
  vector<int> added(m.addCount() + 1);
  copy(m.getAdded().begin(), m.getAdded().end(), added.begin());
  added.back() = getAttrId();
  retval.push_back(added);
  return retval;
}

list<vector<int> > ChangeAttribute::makeAttr
(Move& m, Solution& solution, volumeType, penaltyType) const
{
  list<vector<int> > retval;
  if(m.addCount() != m.removeCount())
    errorExit("ChangeAttribute::makeAttr: move not a change");
  Move::c_indexIterator iI = m.getRemoveIndexes().begin();
  Move::c_patchIterator pI = m.getAdded().begin();
  for(int i = 0; i < m.addCount(); ++i, ++iI, ++pI)
    retval.push_back(makeVector
		     (int(solution.getPatch(*iI)),
		      int(*pI), 
		      getAttrId()));
  return retval;
}

list<vector<int> > ChangeAttribute::makeTestAttr
(Move& m, Solution& solution, volumeType, penaltyType) const
{
  list<vector<int> > retval;
  if(m.addCount() != m.removeCount())
    errorExit("ChangeAttribute::makeTestAttr: move not a change");
  Move::c_indexIterator iI = m.getRemoveIndexes().begin();
  Move::c_patchIterator pI = m.getAdded().begin();
  for(int i = 0; i < m.addCount(); ++i, ++iI, ++pI)
    retval.push_back(makeVector
		     (int(*pI),
		      int(solution.getPatch(*iI)),
		      getAttrId()));
  return retval;
}

list<vector<int> > IndexAttribute::makeAttr
(Move& m, Solution&, volumeType, penaltyType) const
{
  list<vector<int> > retval;
  for(Move::c_indexIterator I=m.getRemoveIndexes().begin();
      I != m.getRemoveIndexes().end(); ++I)
    retval.push_back(makeVector(int(*I), getAttrId()));
  return retval;
}

list<vector<int> > PenaltyVolumeLoopAttribute::makeAttr
(Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const
{
  for(int i = 0; i < 2; ++i) lastValues.pop_front();
  lastValues.push_back(int(newVolume));
  lastValues.push_back(int(newPenalty));
  list<vector<int> > retval;
  vector<int> tmp(checkLen*2 + 1);
  copy(lastValues.begin(), lastValues.end(), tmp.begin());
  tmp.back() = getAttrId();
  retval.push_back(tmp);
  return retval;
}

list<vector<int> > PenaltyVolumeLoopAttribute::makeTestAttr
(Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const
{
  vector<int> tmp(checkLen*2 + 1);
  copy(lastValues.begin() + 2, lastValues.end(), tmp.begin());
  tmp[checkLen*2 - 2] = int(newVolume);
  tmp[checkLen*2 - 1] = int(newPenalty);
  tmp[checkLen*2] = getAttrId();
  list<vector<int> > retval;
  retval.push_back(tmp);
  return retval;
}

list<vector<int> > CompositeAttribute::makeAttr
(Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const
{
  vector<vector<int> > subs(attrs.size());
  list<vector<int> > tmplist;
  int len = 0;
  for(unsigned i = 0; i < attrs.size(); ++i) {
    tmplist = attrs[i]->makeAttr(m, sol, newVolume, newPenalty);
    if(tmplist.size() != 1)
      errorExit("CompositeAttribute: only singleton subattributes accepted");
    subs[i] = tmplist.front();
    len += subs[i].size() - 1;
  }
  vector<int> tmpvect(++len);
  int pos = 0;
  for(unsigned i = 0; i < attrs.size(); ++i) {
    copy(subs[i].begin(), subs[i].end()-1, tmpvect.begin() + pos);
    pos += subs[i].size() - 1;
  }
  tmpvect[pos] = getAttrId();
  return list<vector<int> > (1, tmpvect);
}

list<vector<int> > CompositeAttribute::makeTestAttr
(Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const
{
  vector<vector<int> > subs(attrs.size());
  list<vector<int> > tmplist;
  int len = 0;
  for(unsigned i = 0; i < attrs.size(); ++i) {
    tmplist = attrs[i]->makeTestAttr(m, sol, newVolume, newPenalty);
    if(tmplist.size() != 1)
      errorExit("CompositeAttribute: only singleton subattributes accepted");
    subs[i] = tmplist.front();
    len += subs[i].size() - 1;
  }
  vector<int> tmpvect(++len);
  int pos = 0;
  for(unsigned i = 0; i < attrs.size(); ++i) {
    copy(subs[i].begin(), subs[i].end()-1, tmpvect.begin() + pos);
    pos += subs[i].size() - 1;
  }
  tmpvect[pos] = getAttrId();
  return list<vector<int> > (1, tmpvect);
}

list<vector<int> > NeighborhoodAttribute::makeAttr
(Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const
{
  return neighborhood.makeAttr(m, sol, newVolume, newPenalty);
}

list<vector<int> > NeighborhoodAttribute::makeTestAttr
(Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const
{
  return neighborhood.makeTestAttr(m, sol, newVolume, newPenalty);
}
