#ifndef RUNLIST_H
#define RUNLIST_H
/*
Overarching structure over several runs.
Keeps track of voxels, speeks w/ RunAction & ParticleGenerator.
Calls RunManager->beamOn(XXX)
Commanded directly from main (after user inputs); is its own messenger.
*/

#include "MapRunAction.hh"
#include "L200ParticleGenerator.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"

class G4VAnalysisManager;

class RunList : public G4UImessenger{
public:
	RunList(L200ParticleGenerator* generator, MapRunAction* mra);	//default ctor
	~RunList();	//for now not virtual; do we need to derive from here? (only sober: don't drink and derive!!!)

	//so far I don't have the messenger used here since voxel stuff runs locally inside the generator
	virtual void SetNewValue(G4UIcommand *cmd, G4String newValue);	//@override G4UImessenger

	void startRuns();	//starts all runs until all voxels in the generator are run through

private:
	L200ParticleGenerator* generator;
	MapRunAction* mra;
	
	G4UIdirectory* writeDir;
  	G4UIcmdWithAString* writeFilename;

	G4String filename;

	void openFile();
	void clearVars();
	void writeRun();	//writes single run to file using the ana manager

	G4int count;
	G4double voxelX;		//voxel middle point
	G4double voxelY;
	G4double voxelZ;

	G4VAnalysisManager* analysis;		//for writing out counts
};

#endif
