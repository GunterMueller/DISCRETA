// solution.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <algo.h>
#include <fstream.h>
#include <strstream.h>

#include "packcover.H"
#include "solution.H"

ostream& operator<<(ostream& out, const Solution& s)
{
  s.print(out);
  return out;
}

void Solution::removePatches(const vector<int>& indexes)
{
  if(indexes.empty()) return;
  vector<int> ix(indexes);
  sort(ix.begin(), ix.end());
  for(vector<int>::reverse_iterator I = ix.rbegin(); I != ix.rend(); ++I)
    removePatch(*I);
}

void Solution::removePatches(const list<int>& indexes)
{
  if(indexes.empty()) return;
  vector<int> ix(indexes.size());
  copy(indexes.begin(), indexes.end(), ix.begin());
  removePatches(ix);
}

void Solution::addPatches(const list<patchType>& patches)
{
  for(list<patchType>::const_iterator I = patches.begin();
      I != patches.end(); ++I)
    addPatch(*I);
}

void Solution::changePatches(const Move& move)
{
  Move::c_patchIterator pI(move.getAdded().begin());
  Move::c_indexIterator iI(move.getRemoveIndexes().begin());
  for(int i = 0; i < move.addCount(); ++i, ++iI, ++pI)
    patches[*iI] = *pI;
}

void Solution::writeToFile(const string fileName)
{
  int fd = creat(fileName.c_str(), S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if(fd == -1) errorExit(string("Solution::writeToFile: cannot open file ") +
			 fileName);
  patchType patch;
  for(int i = 0; i < numPatches(); ++i) {
    patch = getPatch(i);
    if(write(fd, &patch, sizeof(patch)) != sizeof(patch))
      errorExit(string("Solution::writeToFile: write failed ") + fileName);
  }
  if(close(fd) == -1) 
    errorExit(string("Solution::writeToFile: close failed: ") +
	      fileName);
}

void Solution::readFromFile(const string fileName)
{
  Solution::clear();
  int fd = open(fileName.c_str(), O_RDONLY);
  if(fd == -1) errorExit(string("Solution::readFromFile: cannot open file ") +
			 fileName);
  patchType patch;
  ssize_t readBytes;
  while(true) {
    readBytes = read(fd, &patch, sizeof(patch));
    if(readBytes == 0) break;
    if(readBytes != sizeof(patch))
      errorExit(string("Solution::readFromFile: read failed ") + fileName);
    Solution::addPatch(patch);
  }
  if(numPatches() == 0)
    errorExit(string("Solution::readFromFile: nothing to read: ") +
	      fileName);
  if(close(fd) == -1) 
    errorExit(string("Solution::readFromFile: close failed: ") +
	      fileName);
  recalculateThings();
}

void Solution::check()
{
  volumeType vol = getCurrentVolume();
  penaltyType pen = getCurrentPenalty();
  recalculateThings();
  if(vol != getCurrentVolume() || pen != getCurrentPenalty())
    errorExit("volume/penalty check failed");
}

void Solution::printDiscretaOutput
(const string filename, const DataProvider *dp) const
{
  ofstream file(filename.c_str(), ios::out | ios::trunc);
  if(!file) errorExit(string("could not open file ") + filename +
			  string(" for writing"));
  multiset<patchType> ps;
  for(vector<patchType>::const_iterator I = patches.begin();
      I != patches.end(); ++I)
    ps.insert(*I);
  for(patchType p = 0; p < dp->totalPatchCount(); ++p)
    file << char('0' + ps.count(p));
  file << '\n';
}

void CoverSolution::recalculateCoverings()
{
  fill(covered->begin(), covered->end(), lambdaType(0));
  const vector<pointType>* points;
  const vector<lambdaType>* lambdas;
  vector<pointType>::const_iterator pI;
  vector<lambdaType>::const_iterator lI;
  bool sorted;
  for(int i = 0; i < numPatches(); ++i) {
    matrix.coveredBy(getPatch(i), points, lambdas, sorted);
    for(pI = points->begin(), lI = lambdas->begin();
	pI != points->end(); ++pI, ++lI)
      (*covered)[*pI] += *lI;
  }
}

void CoverSolution::recomputePenalty(int extraPatches = 0)
{
  setCurrentPenalty(penaltyType(0));
  if(penalty.isUniform()) {
    penaltyType *penTbl = 
      penalty.penaltyTable(-coverMax, 
			   (numPatches()+extraPatches)*
			   dataProvider.maxCovering());
    vector<lambdaType>::const_iterator cI;
    vector<lambdaType>::const_iterator gI;
    for(cI = (*covered).begin(), gI = (*coverGoal).begin();
	cI != (*covered).end(); ++cI, ++gI)
      addPenalty(*(penTbl + (*cI - *gI)));
  } else {
    int i = 0;
    for(vector<lambdaType>::const_iterator I = (*covered).begin();
	I != (*covered).end(); ++I, ++i)
      addPenalty(penalty.penalty(i, *I, (*coverGoal)[i]));
  }
}

penaltyType CoverSolution::pointPenalty(pointType point) const
{
  if(penalty.isUniform()) {
    penaltyType *penTbl = penalty.penaltyTable(-1,0);
    // assume that sufficiently large table has already been
    // allocated, because pointPenalty refers to the penalty of the
    // current solution, whose penalty has already been computed
    return *(penTbl + (*covered)[point] - (*coverGoal)[point]);
  } else
    return penalty.penalty(point, (*covered)[point], (*coverGoal)[point]);
}

void CoverSolution::recomputeVolume()
{
  currVolume = 0;
  for(int i = 0; i < numPatches(); ++i)
    currVolume += dataProvider.getVolume(getPatch(i));
}

struct MergeRec {
  int ix;
  pointType point;
  bool operator<(const MergeRec& rhs) const { return rhs.point < point; }
  // negated form because we want to make a heap with the smallest
  // element on top
};

void CoverSolution::whatIfChange(Move& move, volumeType& volChange, 
				 penaltyType& penChange)
  // note: the merge loops in this function access a STL vector past
  // the last element. This should be ok in usual STL implementations
  // if the capacity of the vector is large enough (made by reserve or
  // resize calls). However, if the STL implementation performs
  // (non-standard) checks, the sizes of the vectors should be
  // increased by one before the merge loops and decreased by one
  // after them.
{
  matrix.lockCache();
  volChange = 0;
  const int remCount = move.removeCount(), addCount = move.addCount();
  const int changeCount = remCount + addCount;
  vector<bool> sorted(changeCount);
  bool fsorted;
  vector<int> changeCoeff(changeCount);
  typedef const vector<pointType>* pvptr;
  typedef const vector<lambdaType>* lbptr;
  vector<pvptr> points(changeCount);
  vector<lbptr> lambdas(changeCount);
  vector<pointType>::const_iterator ptr[addCount + remCount];
  vector<lambdaType>::const_iterator lptr[addCount + remCount];
  int i = 0;
  for(Move::c_indexIterator I = move.getRemoveIndexes().begin();
      I != move.getRemoveIndexes().end(); ++I, ++i) {
    matrix.coveredBy(getPatch(*I), points[i], lambdas[i], fsorted);
    sorted[i] = fsorted;
    changeCoeff[i] = -1;
    volChange -= dataProvider.getVolume(getPatch(*I));
  }
  for(Move::c_patchIterator I = move.getAdded().begin();
      I != move.getAdded().end(); ++I, ++i) {
    matrix.coveredBy(*I, points[i], lambdas[i], fsorted);
    sorted[i] = fsorted;
    changeCoeff[i] = 1;
    volChange += dataProvider.getVolume(*I);
  }
  int pointSum = 0;
  for(int i = 0; i < changeCount; ++i) {
    ptr[i] = points[i]->begin();
    lptr[i] = lambdas[i]->begin();
    pointSum += points[i]->size();
  }
  penChange = 0;
  penaltyType *penTbl = 
    penalty.isUniform() ? 
    penalty.penaltyTable(-coverMax, 
      (numPatches()+move.addCount()) * dataProvider.maxCovering())
    : 0;
  penaltyType *penPtr;

  if(changeCount == 1) { // single add or remove
    vector<pointType>::const_iterator pI;
    vector<lambdaType>::const_iterator lI;
    if(penalty.isUniform())
      for(pI = ptr[0], lI = lptr[0]; pI != points[0]->end(); 
	  ++pI, ++lI) {
	penPtr = penTbl + (*covered)[*pI] - (*coverGoal)[*pI];
	penChange += *(penPtr + changeCoeff[0]*(*lI)) - *penPtr;
      }
    else
      for(pI = ptr[0], lI = lptr[0]; pI != points[0]->end(); 
	  ++pI, ++lI)
	penChange += 
	  penalty.penalty(*pI, (*covered)[*pI]+changeCoeff[0]*(*lI), 
			  (*coverGoal)[*pI]) -
	  penalty.penalty(*pI, (*covered)[*pI], 
			  (*coverGoal)[*pI]);

  } else { // merge several adds/removes

    if(find(sorted.begin(), sorted.end(), false) != sorted.end()) { 
      // can't merge, do the slow sweep
      for(int i = 0; i < changeCount; ++i)
	updateCovered(points[i], lambdas[i], changeCoeff[i]);
      penaltyType oldPenalty = getCurrentPenalty();
      recomputePenalty(move.addCount());
      penChange = getCurrentPenalty() - oldPenalty;
      for(int i = 0; i < changeCount; ++i)
	updateCovered(points[i], lambdas[i], -changeCoeff[i]);
      setCurrentPenalty(oldPenalty);

    } else  // merge
      if(changeCount == 2) { // merge two
	vector<pointType>::const_iterator I0 = points[0]->begin(), 
	  I1 = points[1]->begin();
	vector<lambdaType>::const_iterator lI0 = lambdas[0]->begin(),
	  lI1 = lambdas[1]->begin();
	if(penalty.isUniform())
	  for(int i = 0; i < pointSum; ++i) 
	    if(*I0 == *I1) {
	      penPtr = penTbl+(*covered)[*I0]-(*coverGoal)[*I0];
	      penChange += *(penPtr + changeCoeff[0]*(*lI0)+
			     changeCoeff[1]*(*lI1)) - *penPtr;
	      ++lI0; ++lI1;
	      ++I0; ++I1;
	      ++i;
	    } else if(*I0 < *I1) {
	      penPtr = penTbl+(*covered)[*I0]-(*coverGoal)[*I0];
	      penChange += *(penPtr + changeCoeff[0]*(*lI0)) - *penPtr;
	      ++lI0;
	      ++I0;
	    } else { // *I1 < *I0
	      penPtr = penTbl+(*covered)[*I1]-(*coverGoal)[*I1];
	      penChange += *(penPtr + changeCoeff[1]*(*lI1)) - *penPtr;
	      ++lI1;
	      ++I1;
	    } 
	else { // non-uniform penalty
	  for(int i = points[0]->size() + points[1]->size(); 
	      i > 0; --i) 
	    if(*I0 == *I1) {
	      penChange += 
		penalty.penalty(*I0, (*covered)[*I0]+changeCoeff[0]*(*lI0)+
				changeCoeff[1]*(*lI1),
				(*coverGoal)[*I0])-
		penalty.penalty(*I0, (*covered)[*I0],
				(*coverGoal)[*I0]);
	      ++lI0; ++lI1;
	      ++I0; ++I1;
	      --i;
	    } else if(*I0 < *I1) {
	      penChange += 
		penalty.penalty(*I0, (*covered)[*I0]+changeCoeff[0]*(*lI0),
				(*coverGoal)[*I0])-
		penalty.penalty(*I0, (*covered)[*I0],
				(*coverGoal)[*I0]);
	      ++lI0;
	      ++I0;
	    } else { // *I1 < *I0
	      penChange += 
		penalty.penalty(*I1, (*covered)[*I1]+changeCoeff[1]*(*lI1),
				(*coverGoal)[*I1])-
		penalty.penalty(*I1, (*covered)[*I1],
				(*coverGoal)[*I1]);
	      ++lI1;
	      ++I1;
	    } 
	}

      } else { // merge many
	vector<MergeRec> mheap(changeCount);
	MergeRec curr, next;
	int mergeCnt, pointChange, ix;
	for(int i = 0; i < changeCount; ++i) {
	  mheap[i].ix = i;
	  mheap[i].point = *ptr[i];
	}
	make_heap(mheap.begin(), mheap.end());
	for(int i = 0; i < pointSum; i += mergeCnt) {
	  mergeCnt = 0;
	  curr = mheap[0];
	  pointChange = 0;
	  while(mheap[0].point == curr.point) {
	    ix = mheap[0].ix;
	    ptr[ix]++; next.ix = ix; next.point = *ptr[ix];
	    pop_heap(mheap.begin(), mheap.end());
	    mheap[addCount+remCount-1] = next;
	    push_heap(mheap.begin(), mheap.end());
	    ++mergeCnt;
	    pointChange += *lptr[ix] * changeCoeff[ix]; 
	    ++lptr[ix];
	  }
	  if(penalty.isUniform()) {
	    penPtr = penTbl + (*covered)[curr.point] - 
	      (*coverGoal)[curr.point];
	    penChange += *(penPtr + pointChange) - *penPtr;
	  } else
	    penChange += 
	      penalty.penalty(curr.point, (*covered)[curr.point] + pointChange,
			      (*coverGoal)[curr.point]) -
	      penalty.penalty(curr.point, (*covered)[curr.point],
			      (*coverGoal)[curr.point]);
	}
      }
  }
  matrix.releaseCache();
}

void CoverSolution::updateCovered
(const vector<pointType>* points, const vector<lambdaType>* lambdas,
 const int addRem)
{
  vector<pointType>::const_iterator poI;
  vector<lambdaType>::const_iterator lI;
  for(poI = (*points).begin(), lI = (*lambdas).begin(); 
      poI != (*points).end(); ++poI, ++lI)
    (*covered)[*poI] += (*lI)*addRem;
}

void CoverSolution::makeMoveVP(Move& move, volumeType volCh, penaltyType penCh)
{
  const vector<pointType>* points;
  const vector<lambdaType>* lambdas;
  bool sorted;
  for(Move::c_indexIterator iI = move.getRemoveIndexes().begin();
      iI != move.getRemoveIndexes().end(); ++iI) {
    matrix.coveredBy(getPatch(*iI), points, lambdas, sorted);
    updateCovered(points, lambdas, -1);
  }
  for(Move::c_patchIterator paI = move.getAdded().begin();
      paI != move.getAdded().end(); ++paI) {
    matrix.coveredBy(*paI, points, lambdas, sorted);
    updateCovered(points, lambdas, 1);
  }    
  Solution::makeMoveVP(move, volCh, penCh);
}

void Solution::initSolAtMostVol(const DataProvider& dp, 
				const volumeType maxvol)
{
  clear();
  volumeType currVol = getCurrentVolume(), vol;
  patchType patch;
  while(true) {
    Move m;
    penaltyType pen;
    patch = dp.getRandomPatch();
    m.addAddPatch(patch);
    whatIfChange(m, vol, pen);
    if((currVol += vol) > maxvol)
      break;
    makeMove(m);
  }
}

void Solution::initSolAtLeastVol(const DataProvider& dp, 
				 const volumeType maxvol)
{
  clear();
  while(true) {
    Move m;
    patchType patch = dp.getRandomPatch();
    m.addAddPatch(patch);
    makeMove(m);
    if(getCurrentVolume() >= maxvol)
      break;
  }
}

void CoverSolution::initSolAtLeastTotalCoverage(int reqCov)
{
  clear();
  int currCov = 0;
  const vector<pointType> *pts;
  const vector<lambdaType> *lbd;
  bool sorted;
  while(currCov < reqCov) {
    Move m;
    patchType patch = dataProvider.getRandomPatch();
    m.addAddPatch(patch);
    makeMove(m);
    matrix.coveredBy(patch, pts, lbd, sorted);
    if(lbd->size() == 1)
      currCov += lbd->front() * pts->size();
    else
      currCov += accumulate(lbd->begin(), lbd->end(), 0);
  }
}

void Solution::initSolAtMostPen
(const DataProvider& dp, const penaltyType maxPen)
  // N.B. Will fill the memory, if the penalty does not increase as
  // the number of patches increases
{
  clear();
  vector<patchType> remaining(dp.totalPatchCount());
  for(int i = 0; i < dp.totalPatchCount(); ++i)
    remaining[i] = patchType(i);
  Move m;
  int sel;
  penaltyType penCh;
  volumeType volCh;
  while(!remaining.empty()) {
    sel = rnd(remaining.size());
    m.clear();
    m.addAddPatch(remaining[sel]);
    whatIfChange(m, volCh, penCh);
    if(getCurrentPenalty() + penCh <= maxPen)
      makeMoveVP(m, volCh, penCh);
    else {
      remaining[sel] = remaining.back();
      remaining.resize(remaining.size() - 1);
    }
  }
}

void Solution::initSolPartitions(const DataProvider& dp,
				 const int partSize)
  // takes one random patch of each partition
{
  clear();
  if(dp.totalPatchCount() % partSize != 0) errorExit("inequal partitions");
  int partCount = dp.totalPatchCount() / partSize;
  for(int i = 0; i < partCount; ++i)
    addPatch(patchType(i*partSize + rnd(partSize)));
  recalculateThings();
}

void Solution::initSolVolumes(const DataProvider& dp, const string arg)
{
  clear();
  istrstream istr(arg.c_str(), arg.length());
  int vol, count;
  while(istr >> vol) {
    if(!(istr >> count))
      errorExit("invalid volume specification for initial solution");
    for(int ci = 0; ci < count; ++ci)
      addPatch(dp.getRandomPatch(vol));
  }
  recalculateThings();
}
