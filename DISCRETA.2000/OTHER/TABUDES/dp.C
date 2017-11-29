// dp.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include "packcover.H"
#include "dp.H"

void DataProvider::computeVolVectors()
{
  map<volumeType,int> counts;
  for(int i = 0; i < totalPatchCount(); ++i)
    counts[volumes[i]]++;
  for(map<volumeType,int>::iterator I = counts.begin(); 
      I != counts.end(); ++I) {
    volVectors[I->first] = new vector<patchType>();
    volVectors[I->first]->reserve(I->second);
  }
  for(int i = 0; i < totalPatchCount(); ++i)
    volVectors[volumes[i]]->push_back(patchType(i));
//    for(hash_map<volumeType,vector<patchType>* >::const_iterator 
//  	I = volVectors.begin(); I != volVectors.end(); ++I)
//      cout << "volume " << I->first << ": " << *I->second << '\n';
}
