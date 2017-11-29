// tenure.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include "packcover.H"
#include "tenure.H"
#include "limits.H"
#include "warning.H"

CompareWindows::CompareWindows(int wl) : windowLength(wl)
{
  if(windowLength < 10)
    warning->print("window lengths smaller than 10 are not recommended,\n"
		   "normal approximation may be significantly inaccurate");
  sample.resize(2);
  for(int i = 0; i < windowLength; ++i)
    for(int j = 0; j < 2; ++j)
      sample[j].push_back(numeric_limits<penaltyType>::max());
  for(int j = 0; j < 2; ++j)
    counts[pair<penaltyType,int>(numeric_limits<penaltyType>::max(), j)] = 
      windowLength;
}

bool CompareWindows::changed()
{
  int rankBase = 1;
  double rankSum = 0.;
  count_iterator I, J;
  I = counts.begin();
  while(I != counts.end()) {
    J = I;
    int clusterCount = 0, cluster0Count = 0;
    while(J != counts.end() && J->first.first == I->first.first) {
      clusterCount += J->second;
      if(J->first.second == 1)
	cluster0Count += J->second;
      ++J;
    }
    rankSum += cluster0Count * (2*rankBase + clusterCount - 1) / 2.;
    rankBase += clusterCount;
    I = J;
  }
  // normal approximation, p = 0.25
  //cout << "rankSum: " << rankSum << '\n';
  double statistic = (rankSum - windowLength*(2.*windowLength+1.)/2.) /
    sqrt(windowLength*windowLength*(2.*windowLength+1.)/12.);
  //cout << "statistic: " << statistic << '\n';
  return fabs(statistic) > 0.67;
}

void CompareWindows::add(penaltyType pen)
{
  penaltyType change = sample[1].front(), remove = sample[0].front();
  sample[1].pop_front();
  sample[1].push_back(pen);
  ++counts[pair<penaltyType,int>(pen,1)];
  sample[0].pop_front();
  sample[0].push_back(change);
  count_iterator I = counts.find(pair<penaltyType,int>(change, 1));
  if(!--I->second)
    counts.erase(I);
  ++counts[pair<penaltyType,int>(change,0)];
  I = counts.find(pair<penaltyType,int>(remove,0));
  if(!--I->second)
    counts.erase(I);
}

void DynamicTenure::newPenalty(penaltyType penalty)
{
  registerPenalty(penalty);
  if(!refreshWindow && stagnated()) {
    if(verbose)
      cout << "stagnated: changing tabu tenure from " << currTenure;
    tenIx = (tenIx + 1) % tenures.size();
    currTenure = tenures[tenIx];
    if(verbose)
      cout << " to " << currTenure << '\n';
    refreshWindow = windowLength + currTenure;
  }
}

void DynamicTenure::registerPenalty(penaltyType penalty)
{
  win->add(penalty);
  if(refreshWindow) --refreshWindow;
}
