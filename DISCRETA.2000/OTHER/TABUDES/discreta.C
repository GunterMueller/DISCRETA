// discreta.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include "packcover.H"
#include "discreta.H"
#include "args.H"
#include "limits.H"
#include <fstream.h>

DiscretaDataProvider::DiscretaDataProvider
(const string infile, const long maxmem)
{
  ifstream from(infile.c_str());
  if(!from)
    errorExit(string("could not read file: ") + infile);
  read(from, maxmem);
  coverGoal = new vector<lambdaType>
    (rows, options->get("lambda", lambdaType(1)));
}

void DiscretaDataProvider::read(istream& istr, long maxmem)
{
  int typeChar;

  // preamble

  int order = -1;
  while(!istr.eof() && istr.good()) {
    typeChar = istr.get();
    if(typeChar == EOF)
      errorExit("unexpected end of discretafile");
    else if(typeChar == '%' || typeChar == '\n') {
      int i;
      if(typeChar == '%') {
	do {
	  typeChar = istr.get();
	} while(typeChar != EOF && isspace(typeChar) && typeChar != '\n');
	if(typeChar != '\n' && typeChar != EOF) {
	  istr.putback(typeChar);
	  string word;
	  istr >> word;
	  if(word == string("order:"))
	    istr >> order;
	} else if(typeChar == '\n')
	  istr.putback(typeChar);
      }
      if(istr.good())
	do { i = istr.get(); } while(istr.good() && i != EOF && i != '\n');
    }
    else { // end of preamble 
      istr.putback(typeChar);
      break;
    }
  }
  if(!istr.good())
    errorExit("Error while reading the preamble.");

  rows = cols = 0;
  istr >> rows >> cols;
  if(!rows || !cols)
    errorExit("could not read the size of the matrix");
  if(rows > numeric_limits<pointType>::max())
    errorExit("point type overflow");
  if(cols > numeric_limits<patchType>::max())
    errorExit("patch type overflow");

  // reads the whole matrix into memory, which isn't very good thing
  // to do with large sparse matrixes

  int elem;
  const int maxLambda = numeric_limits<lambdaType>::max();
  lambdaType *matrix = new lambdaType[rows*cols];
  pointType sentinel = numeric_limits<pointType>::max();
  maxCover = 0;
  for(int i = 0; i < rows; ++i)
    for(int j = 0; j < cols; ++j) {
      istr >> elem;
      if(!istr.good())
	errorExit("unexpected end of discretafile");
      if(elem > maxCover)
	maxCover = elem;
      if(elem > maxLambda || elem < 0)
	errorExit("overflow matrix element");
      matrix[i+j*rows] = lambdaType(elem);
    }
  covers.resize(cols);
  coverTimes.resize(cols);
  for(int j = 0; j < cols; ++j) {
    int elemCount = 0;
    for(int i = 0; i < rows; ++i)
      if(matrix[i+j*rows] > 0)
	elemCount++;
    covers[j].reserve(elemCount+1);
    covers[j].resize(elemCount);
    coverTimes[j].resize(elemCount);
    int ix = 0;
    for(int i = 0; i < rows; ++i)
      if(matrix[i+j*rows] > 0) {
	covers[j][ix] = pointType(i);
	coverTimes[j][ix++] = matrix[i+j*rows];
      }
    covers[j][ix] = sentinel;
  }
  delete[] matrix;
  
  // scan for STABILIZER-ORDER-K-SETS
  string inputrow;
  bool stabsFound = false;
  while(!istr.eof() && istr.good()) {
    istr >> inputrow;
    if(inputrow == string("STABILIZER-ORDER-K-SETS")) {
      stabsFound = true;
      if(order == -1)
	if(!options->getNext("order", order))
	  errorExit("The order of group must be defined in the input file\n"
		    "or by giving option \"order\"");
      break;
    }
  }
  if(stabsFound) {
    for(int i = 0; i < totalPatchCount(); ++i) {
      volumes.resize(totalPatchCount());
      int sord;
      istr >> sord;
      if(!istr.good())
	errorExit("Input error while reading STABILIZER-ORDER-K-SETS");
      volumes[i] = order / sord;
      if(sord * volumes[i] != order)
	errorExit("STABILIZER-ORDER does not divide the order of the group");
    }
  } else {
    warning->print("No stabilizer orders given in file, assuming 1-length orbits");
    volumes = vector<volumeType>(totalPatchCount(), volumeType(1));
  }
  computeVolVectors();
}
