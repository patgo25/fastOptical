#include "globals.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

#include "L200ParticleGenerator.hh"
#include "L200ParticleGeneratorMessenger.hh"

//
//---------------------------------------------------------------------------//

L200ParticleGeneratorMessenger::L200ParticleGeneratorMessenger(L200ParticleGenerator *generator)
: fLiquidArgonGenerator(generator){
	// /MG/generator/LiquidArgon
  fLiquidArgonDirectory = new G4UIdirectory("/generator/");
  fLiquidArgonDirectory->SetGuidance("Set to generate optical photons @128 nm in argon inside cryostat");

  fLiquidArgonSetRadius= new G4UIcmdWithADoubleAndUnit("/generator/SetRadiusMax",this);
  fLiquidArgonSetRadius->SetGuidance("Define Max Radius to Generate points inside a Cylidrical Croystat that has diameter = height");
  fLiquidArgonSetRadius->SetDefaultUnit("cm");
  fLiquidArgonSetRadius->SetUnitCategory("Length");
  fLiquidArgonSetRadius->SetUnitCandidates("micron mm cm m km");

  fLiquidArgonSetRadiusMin= new G4UIcmdWithADoubleAndUnit("/generator/SetRadiusMin",this);
  fLiquidArgonSetRadiusMin->SetGuidance("Define Min Radius to Generate points inside a Cylidrical Croystat that has diameter = height");
  fLiquidArgonSetRadiusMin->SetDefaultUnit("cm");
  fLiquidArgonSetRadiusMin->SetUnitCategory("Length");
  fLiquidArgonSetRadiusMin->SetUnitCandidates("micron mm cm m km");


  fLiquidArgonSetHeight= new G4UIcmdWithADoubleAndUnit("/generator/SetHeight",this);
  fLiquidArgonSetHeight->SetGuidance("Define Max Height to Generate points inside a Cylidrical Croystat that has diameter = height");
  fLiquidArgonSetHeight->SetDefaultUnit("cm");
  fLiquidArgonSetHeight->SetUnitCategory("Length");
  fLiquidArgonSetHeight->SetUnitCandidates("micron mm cm m km");

  fLiquidArgonSetBinWidth= new G4UIcmdWithADoubleAndUnit("/generator/SetBinWidth",this);
  fLiquidArgonSetBinWidth->SetGuidance("Define voxel size for optical maps");
  fLiquidArgonSetBinWidth->SetDefaultUnit("cm");
  fLiquidArgonSetBinWidth->SetUnitCategory("Length");
  fLiquidArgonSetBinWidth->SetUnitCandidates("micron mm cm m km");

  fLiquidArgonSetNParticles = new G4UIcmdWithADouble("/generator/SetNParticles",this);
  fLiquidArgonSetNParticles->SetGuidance("Set number of particles to be generated per event");


  //example
  // /generator//SetCenterVector 0.0 0.0 100.0 cm
  fLiquidArgonSetCenterVector= new G4UIcmdWith3VectorAndUnit("/generator/SetCenterVector",this);
  fLiquidArgonSetCenterVector->SetGuidance("Set Center of generator volume inside a Cylidrical Croystat that has diameter = height");
  fLiquidArgonSetCenterVector->SetGuidance("Default value is (0,0,0) but arrays are not centered on (0,0,0)");
  fLiquidArgonSetCenterVector->SetDefaultUnit("cm");
  fLiquidArgonSetCenterVector->SetUnitCategory("Length");
  fLiquidArgonSetCenterVector->SetUnitCandidates("micron mm cm m km");

  fLiquidArgonSet1D = new G4UIcmdWithABool("/generator/Set1D",this);
  fLiquidArgonSet1D->SetGuidance("Set if the voxel probing is 1 or two dimensional");

  fSetVerboseCmd = new G4UIcmdWithAnInteger("/generator/verbose",this);
  fSetVerboseCmd->SetGuidance("Sets verbosity of generator.");
}


L200ParticleGeneratorMessenger::~L200ParticleGeneratorMessenger()
{
  delete fLiquidArgonDirectory;
  delete fLiquidArgonSetRadius;
  delete fLiquidArgonSetRadiusMin;
  delete fLiquidArgonSetHeight;
  delete fLiquidArgonSetCenterVector;
  delete fLiquidArgonSetBinWidth;
  delete fLiquidArgonSetNParticles;
}

void L200ParticleGeneratorMessenger::SetNewValue(G4UIcommand *cmd, G4String str)
{
  if(cmd == fLiquidArgonSetRadius){
    fLiquidArgonGenerator->SetRadius(fLiquidArgonSetRadius->GetNewDoubleValue(str));
  }
  else if(cmd == fLiquidArgonSetRadiusMin){
    fLiquidArgonGenerator->SetRadiusMin(fLiquidArgonSetRadiusMin->GetNewDoubleValue(str));
  }
  else if(cmd == fLiquidArgonSetHeight){
    fLiquidArgonGenerator->SetHeight(fLiquidArgonSetHeight->GetNewDoubleValue(str));
  }
  else if(cmd == fLiquidArgonSetBinWidth){
    fLiquidArgonGenerator->SetBinWidth(fLiquidArgonSetBinWidth->GetNewDoubleValue(str));
  }
  else if(cmd == fLiquidArgonSetNParticles){
    fLiquidArgonGenerator->SetNParticles(fLiquidArgonSetNParticles->GetNewDoubleValue(str));
  }
  else if(cmd == fLiquidArgonSetCenterVector){
    fLiquidArgonGenerator->SetCenterVector(fLiquidArgonSetCenterVector->GetNew3VectorValue(str));
  }
  else if(cmd == fLiquidArgonSet1D){
    fLiquidArgonGenerator->is1DScan(fLiquidArgonSet1D->GetNewBoolValue(str));
  }else if(cmd == fSetVerboseCmd){
		fLiquidArgonGenerator->setVerbosity(fSetVerboseCmd->GetNewIntValue(str));
	  }
}
