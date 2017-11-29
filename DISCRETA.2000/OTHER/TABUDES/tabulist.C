// tabulist.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include "packcover.H"
#include "tabulist.H"
#include "utils.H"
#include "attr.H"
#include "tenure.H"

void TabuList::purgeTime(int t)
{
  TabuEntry timeCmp(makeVector(0), t);

  typedef multiset<TabuEntry, TabuEntry::TimeCmp>::iterator timeIt;
  pair<timeIt,timeIt> timeMatch = timeSet.equal_range(timeCmp);
  timeIt timeI;

  for(timeI = timeMatch.first; timeI != timeMatch.second; ++timeI) {
    typedef multiset<TabuEntry, TabuEntry::AttribCmp>::iterator attrIt;
    pair<attrIt,attrIt> attrMatch = attrSet.equal_range(*timeI);
    multiset<TabuEntry, TabuEntry::AttribCmp>::iterator remI, remTmp;
    
    for(remI = attrMatch.first; remI != attrMatch.second; ) {
      remTmp = remI++;
      if(remTmp->getRemovalTime() == t)
	attrSet.erase(remTmp);
    }
    // can't delete all entries with these attributes, there might be
    // such entries with larger removalTime
  }

  timeSet.erase(timeMatch.first, timeMatch.second);
};

ostream& operator<<(ostream& out, const TabuEntry& te)
{
  out << "TabuEntry: removalTime = " << te.getRemovalTime() <<
    ", attr = ";
  const vector<int>& v(te.getAttr());
  copy(v.begin(), v.end(), ostream_iterator<int> (out, ","));
  return out;
}

ostream& operator<<(ostream& out, TabuList& tl)
{
//   out << "TabuList (timeSet):\n";
//   copy(tl.timeSet.begin(), tl.timeSet.end(), 
//        ostream_iterator<TabuEntry> (out, "\n"));
  out << "TabuList (attrSet):\n";
  copy(tl.attrSet.begin(), tl.attrSet.end(), 
       ostream_iterator<TabuEntry> (out, "\n"));
  return out;
}

void TabuList::makeMove
(Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty, int now)
{
  purgeOld(now);
  for(list<pair<Attribute *const, Tenure *const> >::iterator 
	attrIter = attributes.begin();
      attrIter != attributes.end(); ++attrIter) {
    attrIter->second->newPenalty(newPenalty);
    list<vector<int> > attrList = 
      (attrIter->first)->makeAttr(m, sol, newVolume, newPenalty);
    for(list<vector<int> >::iterator vecIter = attrList.begin();
	vecIter != attrList.end(); ++vecIter) {
      TabuEntry te(*vecIter, now+attrIter->second->time());
      addTabu(te);
    }
  }
}

void TabuList::skipMove(int now)
{
  purgeOld(now);
}

bool TabuList::isTabuMove
(Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const
{
  for(list<pair<Attribute *const, Tenure *const> >::const_iterator 
	attrIter = attributes.begin();
      attrIter != attributes.end(); ++attrIter) {
    list<vector<int> > attrList = 
      (attrIter->first)->makeTestAttr(m, sol, newVolume, newPenalty);
    for(list<vector<int> >::iterator vecIter = attrList.begin();
	vecIter != attrList.end(); ++vecIter) {
      TabuEntry te(*vecIter, 0);
      if(isTabu(te))
	return true;
    }
  }
  return false;
}
