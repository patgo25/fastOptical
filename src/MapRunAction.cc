

#include "G4Timer.hh"
#include "MapRunAction.hh"
#include "G4Run.hh"
#include "globals.hh"


MapRunAction::MapRunAction(size_t nrOfVolumeIndices)
	: G4UserRunAction(), hitCount(nrOfVolumeIndices)
{

}

MapRunAction::~MapRunAction(){}

void MapRunAction::BeginOfRunAction(const G4Run*){
	for(size_t i = 0; i < hitCount.size(); i++){
		hitCount[i] = 0;	//reset counter
	}
}

void MapRunAction::EndOfRunAction(const G4Run*){
	for(size_t i = 0; i < hitCount.size(); i++){//volID starts with one
		G4cout << "Vol "<<i+1<<" --> "<<hitCount[i] << " counts."<<std::endl;
	}
}

void MapRunAction::increment(G4int volID){
	//can be removed for speedup when we are confident, that no shit is going on
	if(volID <= 0){
		G4cout << "ERROR: volID out of bounds: "<<volID<<G4endl;
		G4Exception("MapRunAction::increment","volIDOutOfBounds",RunMustBeAborted,"volume ID of sensitive volume out of bounds: check /g4simple/setVolID in macro");
	}
	G4int index = volID-1;
	if(index >= hitCount.size()) hitCount.resize(index+1, 0);	//fill up missing intermediates with 0
	hitCount[index]++;
}
