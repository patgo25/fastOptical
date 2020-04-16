#ifndef MARRUNACTION_H
#define MAPRUNACTION_H

#include "globals.hh"
#include "G4UserRunAction.hh"

class G4Timer;
class G4Run;

class MapRunAction : public G4UserRunAction
{
  public:
    MapRunAction(size_t nrOfVolumeIndices);
    virtual ~MapRunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

	void increment(G4int volID);	//increments one hit for a specific volume ID

  private:
    //G4Timer* fTimer;
	std::vector<G4int> hitCount;	//hit count per volume (ID = index in vector)
									//I use this since I think its a bit faster than using a map.
									//however, we should know the size of the vector beforehand.
};













#endif
