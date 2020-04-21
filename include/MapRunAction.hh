#ifndef MARRUNACTION_H
#define MARRUNACTION_H

#include "globals.hh"
#include "G4UserRunAction.hh"

class G4Timer;
class G4Run;

class MapRunAction : public G4UserRunAction
{
  public:
    MapRunAction(size_t nrOfVolumeIndices);//nrOfVolumeIndices not needed as vector gets resized on demand
    virtual ~MapRunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

	void increment(G4int volID);	//increments one hit for a specific volume ID
							//volID convention: >0; (0 not allowed, since that means no ID given). 
							//also < 0 not allowed

	G4int getCount(G4int volID){return hitCount.at(volID-1);};
	G4int getVolumeNr() {return hitCount.size();};

  private:
    //G4Timer* fTimer;
	std::vector<G4int> hitCount;	//hit count per volume ((volID-1) = index in vector)
									//I use this since I think its a bit faster than using a map.
									//however, we should know the size of the vector beforehand.
};













#endif
