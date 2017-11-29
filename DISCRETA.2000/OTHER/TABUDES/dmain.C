// dmain.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

// main program for tabudes

#include <iostream.h>
#include <list.h>

#include "packcover.H"
#include "main.H"

#include "args.H"
#include "checked.H"
#include "limits.H"

#include "solution.H"
#include "tabulist.H"
#include "neighborhood.H"
#include "tenure.H"
#include "algorithms.H"

#include "discreta.H"

#include "warning.H"

OptionLine *options = 0;
int verbose;
Warning *warning = 0;

int main(int argc, char **argv)
{
  OptionLine opt(argc, argv);
  options = &opt;
  Warning warn(cout);
  warning = &warn;
  bool solutionFound = false;
  
  long rseedval;
  if(opt.getNext("seed", rseedval))
    setSeed(rseedval);
  else {
    int seed = time(NULL);
    cout << "PRNG seed: " << seed << '\n';
    srandom(seed);
  }

  list<Neighborhood*> neighborhoods;
  vector<Neighborhood*> neighborhood;
  CoverPenalty *coverPenalty = 0;
  enum { cover, pack } probType = cover;

  verbose = opt.get("verbose", 2);
  string prob("discreta");
  probType = cover;

  CoveringDataProvider *cdp = 0;
  CoverSolution *csol = 0;
  CoverMatrix *cmat = 0;
  cdp = new DiscretaDataProvider
    (opt.getMust("discretafile", string()), 
     opt.get("maxmem", numeric_limits<int>::max()));
  string covPenDflt("design");
  if(opt.get("coverpenalty", covPenDflt) == string("design"))
    coverPenalty = new DesignPenalty();
  else if(opt.get("coverpenalty", covPenDflt) == string("asymhigh"))
    coverPenalty = new AsymHighDesignPenalty
      (opt.getMust("asymadd", penaltyType()));
  else if(opt.get("coverpenalty", covPenDflt) == string("asymlow"))
    coverPenalty = new AsymLowDesignPenalty
      (opt.getMust("asymadd", penaltyType()));
  else if(opt.get("coverpenalty", covPenDflt) == string("cover"))
    coverPenalty = new CoveringPenalty();
  else if(opt.get("coverpenalty", covPenDflt) == string("pack"))
    coverPenalty = new PackingPenalty();
  else errorExit(string("Unknown penalty type: ") + 
		 opt.getMust("coverpenalty", string()));
  cmat = new CoverMemoryMatrix(*cdp, opt.get
			       ("maxmem", numeric_limits<int>::max()));
  csol = new CoverSolution(*cmat, *cdp, *coverPenalty);
  
  // initial solution
  string initsol = opt.get("initsol", string("coverage"));
  if(initsol == string("atleastvol"))
    csol->initSolAtLeastVol(*cdp, opt.getMust("vol", int()));
  else if(initsol == string("atmostvol"))
    csol->initSolAtMostVol(*cdp, opt.getMust("vol", int()));
  else if(initsol == string("read"))
    csol->readFromFile(opt.getMust("infile", string()));
  else if(initsol == string("atmostpen"))
    csol->initSolAtMostPen(*cdp, 
			       opt.getMust("pen", penaltyType(0)));
  else if(initsol == string("coverage"))
    csol->initSolAtLeastTotalCoverage
      (opt.getMust("lambda", int()) * cdp->totalPointCount());
  else if(initsol == string("volumes"))
    csol->initSolVolumes(*cdp, opt.getMust("volcounts", string()));
  else if(initsol == string("empty"))
      ; // nothing needs to be done
  else
    errorExit(string("unknown initial solution type: ") + initsol);

  // neighborhood

  PopTopStack<NeighParseStackElem> neighParseStack;

  string neighParseStr = opt.get("neigh", string());
  istrstream istr(neighParseStr.c_str(), neighParseStr.length());
  string token;
  while(!istr.eof() && neighParseStr != string()) {
    istr >> token;
    Neighborhood *neigh = 0;
    if(token == string("changeone"))
      neigh = new ChangeManyNeighborhood(*csol, *cdp, 1);
    else if(token == string("vchangeone"))
      neigh = new VChangeOneNeighborhood(*csol, *cdp);
    else if(token == string("addmany"))
      neigh = new AddManyNeighborhood
	(*cdp, atoi(neighParseStack.popTop().getStr().c_str()));
    else if(token == string("addmanyunique"))
      neigh = new AddManyUniqueNeighborhood
	(*cdp, atoi(neighParseStack.popTop().getStr().c_str()));
    else if(token == string("addone"))
      neigh = new AddManyNeighborhood(*cdp, 1);
    else if(token == string("removemany"))
      neigh = new RemoveManyNeighborhood
	(*csol, atoi(neighParseStack.popTop().getStr().c_str()));
    else if(token == string("removeone"))
      neigh = new RemoveManyNeighborhood(*csol, 1);
    else if(token == string("changemany"))
      neigh = new ChangeManyNeighborhood
	(*csol, *cdp, 
	 atoi(neighParseStack.popTop().getStr().c_str()));
    else if(token == string("changemanyunique"))
      neigh = new ChangeManyUniqueNeighborhood
	(*csol, *cdp, 
	 atoi(neighParseStack.popTop().getStr().c_str()));
    else if(token == string("empty"))
      neigh = new EmptyNeighborhood();
    else if(token == string("compose")) {
      int n = atoi(neighParseStack.popTop().getStr().c_str());
      CompositeNeighborhood *compNeigh;
      neigh = compNeigh = new CompositeNeighborhood();
      for(int i = 0; i < n; ++i)
	compNeigh->addNeighborhood(neighParseStack.popTop().getNptr());
    }
    else if(token == string("union")) {
      int n = atoi(neighParseStack.popTop().getStr().c_str());
      AdditiveNeighborhood *addNeigh;
      neigh = addNeigh = new AdditiveNeighborhood();
      for(int i = 0; i < n; ++i)
	addNeigh->addNeighborhood(neighParseStack.popTop().getNptr());
    }
    else if(token == string("correct")) {
      Neighborhood *incN, *decN;
      decN = neighParseStack.popTop().getNptr();
      incN = neighParseStack.popTop().getNptr();
      neigh = new CorrectCoveringNeighborhood
	(*csol, incN, decN, *cmat);
    }
    else if(token == string("correctrandom")) {
      Neighborhood *incN, *decN;
      decN = neighParseStack.popTop().getNptr();
      incN = neighParseStack.popTop().getNptr();
      neigh = new CorrectRandomCoveringNeighborhood
	(*csol, incN, decN, *cmat);
    }
    else if(token == string("largecorrect")) {
      Neighborhood *incN, *decN;
      decN = neighParseStack.popTop().getNptr();
      incN = neighParseStack.popTop().getNptr();
      neigh = new LargeCorrectCoveringNeighborhood
	(*csol, incN, decN, *cmat);
    }
    else if(token == string("balance")) {
      Neighborhood *incN, *decN;
      decN = neighParseStack.popTop().getNptr();
      incN = neighParseStack.popTop().getNptr();
      neigh = new BalanceCoveringNeighborhood(*csol, incN, decN, *cmat);
    }
    else if(token == string("largebalance")) {
      Neighborhood *incN, *decN;
      decN = neighParseStack.popTop().getNptr();
      incN = neighParseStack.popTop().getNptr();
      neigh = new LargeBalanceCoveringNeighborhood(*csol, incN, decN, *cmat);
    }
    else if(token == string("select")) {
      int goal = atoi(neighParseStack.popTop().getStr().c_str());
      int n = atoi(neighParseStack.popTop().getStr().c_str());
      MultiNeighborhood *multi;
      neigh = multi = new MultiNeighborhood(goal);
      for(int i = 0; i < n; ++i)
	multi->addNeighborhood(neighParseStack.popTop().getNptr());
    }
    else
      neighParseStack.push(NeighParseStackElem(token));
    if(neigh) {
      neighborhoods.push_back(neigh);
      neighParseStack.push(NeighParseStackElem(neigh));
    }
  }
  while(!neighParseStack.empty())
    neighborhood.push_back(neighParseStack.popTop().getNptr());
  reverse(neighborhood.begin(), neighborhood.end());

  string algo = opt.get("algo", string());
  if(neighborhood.empty() && algo != string("none")) 
    errorExit("unspecified neighborhood");

  // attributes & tabu list?
  vector<Attribute*> attributes;
  TabuList *tabuList = 0;
  vector<Tenure*> tenures;
  if(algo == string("tabu")) {
    tabuList = new TabuList();
    PopTopStack<AttrParseStackElem> attrParseStack;
    string attrParseStr = opt.getMust("attr", string());
    istrstream istr(attrParseStr.c_str(), attrParseStr.length());
    string token;
    while(!istr.eof()) {
      istr >> token;
      Attribute *attr = 0;
      if(token == string("remove"))
	attr = new RemoveAttribute();
      else if(token == string("add"))
	attr = new AddAttribute();
      else if(token == string("addall"))
	attr = new AddAllAttribute();
      else if(token == string("index"))
	attr = new IndexAttribute();
      else if(token == string("change"))
	attr = new ChangeAttribute();
      else if(token == string("penvol"))
	attr = new PenaltyVolumeLoopAttribute
	  (atoi(attrParseStack.popTop().getStr().c_str()));
      else if(token == string("and")) {
	int count = atoi(attrParseStack.popTop().getStr().c_str());
	CompositeAttribute* cattr;
	attr = cattr = new CompositeAttribute();
	for(int i = 0; i < count; ++i)
	  cattr->addAttribute(attrParseStack.popTop().getAptr());
      }
      else if(token == string("tl")) {
	int tl = atoi(attrParseStack.popTop().getStr().c_str());
	Tenure *ten = new ConstantTenure(tl);
	tenures.push_back(ten);
	tabuList->addAttributeToList(attrParseStack.popTop().getAptr(), ten);
      }
      else if(token == string("tls")) {
	int count = atoi(attrParseStack.popTop().getStr().c_str());
	int tl = atoi(attrParseStack.popTop().getStr().c_str());
	Tenure *ten = new ConstantTenure(tl);
	for(int i = 0; i < count; ++i)
	  tabuList->addAttributeToList(attrParseStack.popTop().getAptr(), ten);
      }
      else if(token == string("dtl")) {
	vector<int> tenureTimes;
	while(true) {
	  string ten = attrParseStack.popTop().getStr();
	  if(ten == string("dt")) break;
	  tenureTimes.push_back(atoi(ten.c_str()));
	}
	int window = atoi(attrParseStack.popTop().getStr().c_str());
	if(tenureTimes.empty()) errorExit("no tenures in dtl definition");
	reverse(tenureTimes.begin(), tenureTimes.end());
	Tenure *ten = new DynamicTenure(tenureTimes, window);
	tenures.push_back(ten);
	tabuList->addAttributeToList(attrParseStack.popTop().getAptr(), ten);
      }
      else if(token == string("ctl")) {
	vector<int> tenureTimes;
	while(true) {
	  string ten = attrParseStack.popTop().getStr();
	  if(ten == string("ct")) break;
	  tenureTimes.push_back(atoi(ten.c_str()));
	}
	int changeCount = atoi(attrParseStack.popTop().getStr().c_str());
	if(tenureTimes.empty()) errorExit("no tenures in ctl definition");
	reverse(tenureTimes.begin(), tenureTimes.end());
	Tenure *ten = new ShiftTenure(tenureTimes, changeCount);
	tenures.push_back(ten);
	tabuList->addAttributeToList(attrParseStack.popTop().getAptr(), ten);
      }
      else
	attrParseStack.push(AttrParseStackElem(token));
      if(attr) {
	attributes.push_back(attr);
	attrParseStack.push(AttrParseStackElem(attr));
      }
    }
    if(!attrParseStack.empty())
      errorExit("attributes: parse error: final stack not empty");
    if(attributes.empty())
      warning->print("no attributes defined");
  }

  // algorithm
  if(algo == string("SA"))
    basicSimulatedAnnealing(*csol, neighborhood, 
			    opt.get("endlimit", penaltyType(0)));
  else if(algo == string("iterpen"))
    iterateWithinPenaltyLimit
      (*csol, neighborhood,
       opt.get("minpen", numeric_limits<penaltyType>::min()),
       opt.get("maxpen", numeric_limits<penaltyType>::max()),
       opt.get("maxiter", numeric_limits<int>::max()));
  else if(algo == string("tabu"))
    basicTabuSearch(*csol, neighborhood, *tabuList,
		    opt.get("endlimit", penaltyType(0)));
  else if(algo == string("twoway"))
    twoWaySearch(*csol, neighborhood, opt.get("endlimit", penaltyType(0)));
  else if(algo == string("rwalk"))
    randomWalk(*csol, neighborhood, opt.getMust("maxiter", int()));
  else if(algo == string("rwalkpen"))
    randomWalkWithinPenaltyLimit
      (*csol, neighborhood,
       opt.get("minpen", numeric_limits<penaltyType>::min()),
       opt.get("maxpen", numeric_limits<penaltyType>::max()),
       opt.get("maxiter", numeric_limits<int>::max()));
  else if(algo == string("threshold"))
    thresholdAccepting(*csol, neighborhood, 
		       opt.get("endlimit", penaltyType(0)));
  else if(algo == string("record"))
    recordToRecordTravel(*csol, neighborhood, 
			 opt.get("endlimit", penaltyType(0)));
  else if(algo == string("deluge"))
    greatDeluge(*csol, neighborhood, 
		opt.get("endlimit", penaltyType(0)));
  else if(algo == string("local"))
    localOptimize(*csol, neighborhood);
  else if(algo == string("none"))
      ; // no search
  else
    errorExit("no algorithm specified");

  // save eventual results
  if(csol->getCurrentPenalty() <= 
     opt.get("maxpen", numeric_limits<penaltyType>::max()) &&
     csol->getCurrentPenalty() >= 
     opt.get("minpen", numeric_limits<penaltyType>::min()) &&
     csol->getCurrentVolume() <= 
     opt.get("maxvol", numeric_limits<volumeType>::max()) &&
     csol->getCurrentVolume() >= 
     opt.get("minvol", numeric_limits<volumeType>::min())) {
    solutionFound = true;
    if(verbose) {
      cout << "solution found (penalty = " << csol->getCurrentPenalty()
	   << ", volume = " << csol->getCurrentVolume() << "):\n"
	   << *csol << '\n';
    }
    if(opt.get("check", true)) {
      if(verbose) cout << "checking solution...\n";
      csol->check();
      if(verbose) cout << "check passed\n";
    }
    if(opt.get("outfile", string()) != string())
      csol->writeToFile(opt.getMust("outfile", string()));
    if(csol->getCurrentPenalty() <= 
       opt.get("discretapenalty", penaltyType(0)) &&
       opt.get("discretaoutput", string()) != string())
      csol->printDiscretaOutput(opt.getMust("discretaoutput", string()),
				cdp);
  }

  // cleanup
  delete csol;
  for(list<Neighborhood*>::iterator I = neighborhoods.begin();
      I != neighborhoods.end(); ++I)
    delete *I;
  for(vector<Attribute*>::iterator I = attributes.begin();
      I != attributes.end(); ++I)
    delete *I;
  for(vector<Tenure*>::iterator I = tenures.begin();
      I != tenures.end(); ++I)
    delete *I;
  delete cdp;
  delete cmat;
  delete coverPenalty;
  delete tabuList;

  warn.flush();

  return !solutionFound;
}
