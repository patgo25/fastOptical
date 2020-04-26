#ifndef L200FiberPhysics_h
#define L200FiberPhysics_h

#include "G4VPhysicsConstructor.hh"
#include "L200OpBoundaryProcess.hh"

class L200FiberPhysics : public G4VPhysicsConstructor {

public:
	L200FiberPhysics(G4int verbose=0, const G4String& name="L200Op");
	virtual ~L200FiberPhysics();
	virtual void ConstructParticle();
	virtual void ConstructProcess();

	//Magic setter
	void setFiberHitProb(G4double value){theProb = value;}
	void setMagicMaterialName(G4String value){theTPBMagicMaterialName = value;}
	void setLArWL(G4double value){theLArWL = value;}

private:
	L200OpBoundaryProcess* fL200OpBoundaryProcess;
	G4double theProb;
	G4String theTPBMagicMaterialName;
	G4double theLArWL;
};
#endif
