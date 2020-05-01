#ifndef L200DetectorMessenger_h
#define L200DetectorMessenger_h

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"


class L200DetectorConstruction; //Forward declration

class L200DetectorMessenger : public G4UImessenger{

public:
	L200DetectorMessenger(L200DetectorConstruction* thedet);
	virtual ~L200DetectorMessenger();
	virtual void SetNewValue(G4UIcommand* command, G4String value);

private:
	L200DetectorConstruction* det;

	G4UIdirectory* geomDir;
	G4UIdirectory* geomInnerShroudDir;
	G4UIdirectory* geomOuterShroudDir;
	G4UIdirectory* geomWLSRDir;
	G4UIdirectory* geomGeDir;

	G4UIdirectory* opticsDir;

	G4UIcmdWithADoubleAndUnit* innerShroudInnerRadiusCmd;
	G4UIcmdWithADoubleAndUnit* innerShroudOuterRadiusCmd;
	G4UIcmdWithADoubleAndUnit* innerShroudHeightCmd;
	G4UIcmdWithADoubleAndUnit* innerShroudZOffsetCmd;

	G4UIcmdWithADoubleAndUnit* outerShroudInnerRadiusCmd;
	G4UIcmdWithADoubleAndUnit* outerShroudOuterRadiusCmd;
	G4UIcmdWithADoubleAndUnit* outerShroudHeightCmd;
	G4UIcmdWithADoubleAndUnit* outerShroudZOffsetCmd;

	G4UIcmdWithADoubleAndUnit* wlsrRadiusCmd;
	G4UIcmdWithADoubleAndUnit* wlsrHeightCmd;
	G4UIcmdWithADoubleAndUnit* wlsrTPBThicknessCmd;
	G4UIcmdWithADoubleAndUnit* wlsrCuThicknessCmd;
	G4UIcmdWithADoubleAndUnit* wlsrTetraTexThicknessCmd;

	G4UIcmdWithADoubleAndUnit* cryostatWallThicknessCmd;

	G4UIcmdWithADoubleAndUnit* setGeDiscHeightCmd;
	G4UIcmdWithADoubleAndUnit* setGeDiscRadCmd;
	G4UIcmdWithADoubleAndUnit* setGeDiscGapCmd;
	G4UIcmdWithADoubleAndUnit* setGeArrayRadCmd;
	G4UIcmdWithAnInteger* setNrGeDetPerStringCmd;
	G4UIcmdWithAnInteger* setNrGeStringsCmd;

	G4UIcmdWithADoubleAndUnit* lArAbsLengthCmd;
	G4UIcmdWithADoubleAndUnit* visAbsLengthCmd;
	G4UIcmdWithADoubleAndUnit* lArScintWLCmd;
	G4UIcmdWithADoubleAndUnit* tpbScintWLCmd;
	G4UIcmdWithABool*	   lArIsRayCmd;
	G4UIcmdWithABool*	   setBlackWLSRCmd;

	G4UIcmdWithoutParameter*   updateCmd;
};
#endif






