#include "L200FiberPhysics.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"

L200FiberPhysics::L200FiberPhysics(G4int verbose, const G4String& name) : G4VPhysicsConstructor(name) {

	fL200OpBoundaryProcess = NULL;
	verboseLevel = verbose;
	theProb = 0.4;
	theTPBMagicMaterialName="LiquidArgonFiber";
	theLArWL = 128*nm;
}

L200FiberPhysics::~L200FiberPhysics(){

	delete fL200OpBoundaryProcess;
}

void L200FiberPhysics::ConstructParticle() {

	G4OpticalPhoton::OpticalPhotonDefinition();
}

void L200FiberPhysics::ConstructProcess() {

	fL200OpBoundaryProcess = new L200OpBoundaryProcess();

	fL200OpBoundaryProcess->setFiberHitProb(theProb);
	fL200OpBoundaryProcess->setMagicMaterialName(theTPBMagicMaterialName);
	fL200OpBoundaryProcess->setLArWL(theLArWL);

	G4ProcessManager* pm = 0;
	pm = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

	if(!pm){
		std::ostringstream o;
		o << "Optical Photon w/o a Process Manager!";
		G4Exception("L200FiberPhysics::ConstructProcess()","",FatalException, o.str().c_str());
		return;
	}
	fL200OpBoundaryProcess->SetVerboseLevel(verboseLevel);

	pm->AddDiscreteProcess(fL200OpBoundaryProcess);
}
