// matrix.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include <list.h>
#include "packcover.H"
#include "matrix.H"

CoverMemoryMatrix::CoverMemoryMatrix
(CoveringDataProvider& dp, long mxmem) : CoverMatrix(dp, mxmem), memSum(0L)
{
  for(int i = 0; i < dp.totalPatchCount(); ++i)
    if((memSum += dp.memNeeded(patchType(i))) > mxmem)
      errorExit("not enough memory for covering matrix");
  pointVectors.resize(dp.totalPatchCount());
  lambdaVectors.resize(dp.totalPatchCount());
  sortedVector.resize(dp.totalPatchCount());
  memSum = 0L;
  bool fsorted;
  for(int i = 0; i < dp.totalPatchCount(); ++i) {
    dp.getCoverings(patchType(i), pointVectors[i], lambdaVectors[i], fsorted);
    memSum += pointVectors[i]->size() * sizeof(pointType);
    memSum += lambdaVectors[i]->size() * sizeof(lambdaType);
    if(memSum > mxmem)
      errorExit("not enough memory for covering matrix");
    sortedVector[i] = fsorted;
  }
}

void CoverMemoryMatrix::print(ostream& out) const
{
  for(unsigned i = 0; i < pointVectors.size(); ++i)
    out << "patch " << i << " covers: " << *pointVectors[i] << '\n';
}

void CoverMemoryMatrix::calculateInverses() const
{
  if(verbose) cout << "calculating inverse tables...\n";
  typedef list<patchType> patchListType;
  vector<patchListType> coverings(dataProvider.totalPointCount());
  vector<pointType>* points;
  vector<lambdaType>* lambdas;
  bool sorted;
  for(patchType patch = 0; patch < dataProvider.totalPatchCount(); ++patch) {
    coveredBy(patch, points, lambdas, sorted);
    for(pointVectorType::const_iterator I = points->begin();
	I != points->end(); ++I) {
      coverings[*I].push_back(patch);
      memSum += sizeof(patchType);
      if(memSum > getMaxMemory())
	errorExit("out of memory in CoverMemoryMatrix::calculateInverses");
    }
  }
  inverseVectors.resize(dataProvider.totalPointCount());
  for(pointType point = 0; point < dataProvider.totalPointCount(); ++point) {
    inverseVectors[point] = new vector<patchType>(coverings[point].size());
    copy(coverings[point].begin(), coverings[point].end(),
	 inverseVectors[point]->begin());
    coverings[point].clear();
  }
  if(verbose) cout << "inverse tables complete\n";
}
