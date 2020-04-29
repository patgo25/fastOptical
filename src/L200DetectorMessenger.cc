#include "L200DetectorMessenger.hh"
#include "L200DetectorConstruction.hh"

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"

using namespace CLHEP;

L200DetectorMessenger::L200DetectorMessenger(L200DetectorConstruction* thedet) : det(thedet){

	//UI commandsi

	geomDir = new G4UIdirectory("/geometry/");
	geomWLSRDir = new G4UIdirectory("/geometry/wlsr/");
	geomInnerShroudDir = new G4UIdirectory("/geometry/innerShroud/");
	geomOuterShroudDir = new G4UIdirectory("/geometry/outerShroud/");
	opticsDir = new G4UIdirectory("/optics/");

	updateCmd = new G4UIcmdWithoutParameter("/update",this);
	updateCmd->SetGuidance("Update Parameters");
	updateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


	innerShroudInnerRadiusCmd = new G4UIcmdWithADoubleAndUnit("/geometry/innerShroud/innerRadius", this);
	innerShroudInnerRadiusCmd->SetDefaultValue(122.5*mm);
	innerShroudInnerRadiusCmd->SetGuidance("Set the inner radius of the inner fiber shroud");

	innerShroudOuterRadiusCmd = new G4UIcmdWithADoubleAndUnit("/geometry/innerShroud/outerRadius", this);
	innerShroudOuterRadiusCmd->SetDefaultValue(134.5*mm);
	innerShroudOuterRadiusCmd->SetGuidance("Set the outer radius of the inner fiber shroud");

	outerShroudInnerRadiusCmd = new G4UIcmdWithADoubleAndUnit("/geometry/outerShroud/innerRadius", this);
	outerShroudInnerRadiusCmd->SetDefaultValue(283.0*mm);
	innerShroudInnerRadiusCmd->SetGuidance("Set the inner radius of the outer fiber shroud");

	outerShroudOuterRadiusCmd = new G4UIcmdWithADoubleAndUnit("/geometry/outerShroud/outerRadius", this);
	outerShroudOuterRadiusCmd->SetDefaultValue(295.0*mm);
	outerShroudOuterRadiusCmd->SetGuidance("Set the outer radius of the outer fiber shroud");

	wlsrRadiusCmd = new G4UIcmdWithADoubleAndUnit("/geometry/wlsr/radius", this);
	wlsrRadiusCmd->SetDefaultValue(700*mm);
	wlsrRadiusCmd->SetGuidance("Set the radius of the WLSR");

	innerShroudHeightCmd = new G4UIcmdWithADoubleAndUnit("/geometry/innerShroud/height", this);
	innerShroudHeightCmd->SetDefaultValue(1300.0*mm);
	innerShroudHeightCmd->SetGuidance("Set the height of the inner fiber shroud");

	outerShroudHeightCmd = new G4UIcmdWithADoubleAndUnit("/geometry/outerShroud/height", this);
	outerShroudHeightCmd->SetDefaultValue(1500.0*mm);
	outerShroudHeightCmd->SetGuidance("Set the height of the outer fiber shroud");

	wlsrHeightCmd = new G4UIcmdWithADoubleAndUnit("/geometry/wlsr/height", this);
	wlsrHeightCmd->SetDefaultValue(3500.0*mm);
	wlsrHeightCmd->SetGuidance("Set the height of the WLSR");

	wlsrTPBThicknessCmd = new G4UIcmdWithADoubleAndUnit("/geometry/wlsr/tpbThickness", this);
	wlsrTPBThicknessCmd->SetDefaultValue(0.001*mm);
	wlsrTPBThicknessCmd->SetGuidance("Set the TPB coating thickness of the WLSR");

	wlsrCuThicknessCmd = new G4UIcmdWithADoubleAndUnit("/geometry/wlsr/cuThickness", this);
	wlsrCuThicknessCmd->SetDefaultValue(0.03*mm);
	wlsrCuThicknessCmd->SetGuidance("Set Thickness of the copper foil");

	wlsrTetraTexThicknessCmd = new G4UIcmdWithADoubleAndUnit("/geometry/wlsr/tetraTexThickness", this);
	wlsrTetraTexThicknessCmd->SetDefaultValue(0.01*mm);
	wlsrTetraTexThicknessCmd->SetGuidance("Set the thickness of the TetraTex");

	cryostatWallThicknessCmd = new G4UIcmdWithADoubleAndUnit("/geometry/cryostat/wallThickness", this);
	cryostatWallThicknessCmd->SetDefaultValue(10*mm);
	cryostatWallThicknessCmd->SetGuidance("Set the thickness of the cryostat wall");

	lArAbsLengthCmd = new G4UIcmdWithADoubleAndUnit("/optics/lArAbsLength", this);
	lArAbsLengthCmd->SetDefaultValue(20*cm);
	lArAbsLengthCmd->SetGuidance("Set absorption length of LAr scintillation light in LAr");

	visAbsLengthCmd = new G4UIcmdWithADoubleAndUnit("/optics/visAbsLength", this);
	visAbsLengthCmd->SetDefaultValue(1000*m);
	visAbsLengthCmd->SetGuidance("Set absorption length of TPB emission light in LAr");

	lArScintWLCmd = new G4UIcmdWithADoubleAndUnit("/optics/lArScintWL", this);
	lArScintWLCmd->SetDefaultValue(128*nm);
	lArScintWLCmd->SetGuidance("Set LAr scintillation light WL");

	tpbScintWLCmd = new G4UIcmdWithADoubleAndUnit("/optics/tpbScintWL", this);
	tpbScintWLCmd->SetDefaultValue(20*cm);
	tpbScintWLCmd->SetGuidance("Set TPB emission light WL");

	lArIsRayCmd = new G4UIcmdWithABool("/optics/lArRayToggle", this);
	lArIsRayCmd->SetDefaultValue(false);
	lArIsRayCmd->SetGuidance("Set if Ray Scattering can occur in LAr");

	setBlackWLSRCmd = new G4UIcmdWithABool("/optics/setBlackWLSR", this);
	setBlackWLSRCmd->SetDefaultValue(false);
	setBlackWLSRCmd->SetGuidance("true-> make the TPB on the WLSR black");


}

L200DetectorMessenger::~L200DetectorMessenger(){

	
	delete geomInnerShroudDir;
	delete geomOuterShroudDir;
	delete geomWLSRDir;

	delete geomDir;		//delete directory after contents

	delete innerShroudInnerRadiusCmd;
	delete innerShroudOuterRadiusCmd;
	delete outerShroudInnerRadiusCmd;
	delete outerShroudOuterRadiusCmd;
	delete wlsrRadiusCmd;
	delete innerShroudHeightCmd;
	delete outerShroudHeightCmd;
	delete wlsrHeightCmd;
	delete wlsrTPBThicknessCmd;
	delete wlsrCuThicknessCmd;
	delete wlsrTetraTexThicknessCmd;
	delete cryostatWallThicknessCmd;
	delete lArAbsLengthCmd;
	delete visAbsLengthCmd;
	delete lArScintWLCmd;
	delete tpbScintWLCmd;
	delete updateCmd;
	delete lArIsRayCmd;

	delete opticsDir;
}

void L200DetectorMessenger::SetNewValue(G4UIcommand* command, G4String value){
	if(command == innerShroudInnerRadiusCmd){
		det->setinnerShroudInnerR(innerShroudInnerRadiusCmd->GetNewDoubleValue(value));
	}
	if(command == innerShroudOuterRadiusCmd){
		det->setinnerShroudOuterR(innerShroudOuterRadiusCmd->GetNewDoubleValue(value));
	}
	if(command == outerShroudInnerRadiusCmd){
		det->setouterShroudInnerR(outerShroudInnerRadiusCmd->GetNewDoubleValue(value));
	}
	if(command == outerShroudOuterRadiusCmd){
		det->setouterShroudOuterR(outerShroudOuterRadiusCmd->GetNewDoubleValue(value));

	}
	if(command == wlsrRadiusCmd){
		det->setwlsrRadius(wlsrRadiusCmd->GetNewDoubleValue(value));

	}
	if(command == innerShroudHeightCmd){
		det->setinnerShroudHeight(innerShroudHeightCmd->GetNewDoubleValue(value));

	}
	if(command == outerShroudHeightCmd){
	 	det->setouterShroudHeight(outerShroudHeightCmd->GetNewDoubleValue(value));

	}
	if(command == wlsrHeightCmd){
		det->setwlsrHeight(wlsrHeightCmd->GetNewDoubleValue(value));

	}
	if(command == wlsrTPBThicknessCmd){
		det->setwlsrTPBThickness(wlsrTPBThicknessCmd->GetNewDoubleValue(value));

	}
	if(command == wlsrCuThicknessCmd){
		det->setwlsrCuThickness(wlsrCuThicknessCmd->GetNewDoubleValue(value));

	}
	if(command == wlsrTetraTexThicknessCmd){
		det->setwlsrTetraTexThickness(wlsrTetraTexThicknessCmd->GetNewDoubleValue(value));

	}
	if(command == cryostatWallThicknessCmd){
		det->setcryostatWallThickness(cryostatWallThicknessCmd->GetNewDoubleValue(value));

	}
	if(command == lArAbsLengthCmd){
		det->setlArAbsVUV(lArAbsLengthCmd->GetNewDoubleValue(value));

	}

	if(command == visAbsLengthCmd){
		det->setlArAbsVis(visAbsLengthCmd->GetNewDoubleValue(value));

	}
	if(command == lArScintWLCmd){
		det->setlArWL(lArScintWLCmd->GetNewDoubleValue(value));

	}
	if(command == tpbScintWLCmd){
		det->settpbWL(tpbScintWLCmd->GetNewDoubleValue(value));

	}
	if(command == updateCmd){
		det->UpdateGeometry();
	}
	if(command == lArIsRayCmd){
		det->setlArRay(lArIsRayCmd->GetNewBoolValue(value));
	}
	if(command == setBlackWLSRCmd){
		det->setBlackWLSR(setBlackWLSRCmd->GetNewBoolValue(value));
	}
	



}
