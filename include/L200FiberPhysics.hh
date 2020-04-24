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
private:
	L200OpBoundaryProcess* fL200OpBoundaryProcess;
};
#endif
