#ifndef L200DetectorConstruction_h
#define L200DetectorConstruction_h

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "globals.hh"
#include "Randomize.hh"




class G4VPhysicalVolume;
class G4LogicalVolume;


class L200DetectorConstruction : public G4VUserDetectorConstruction{

public:
    	L200DetectorConstruction();
    	virtual ~L200DetectorConstruction();

	//wird vom RunManager aufgerufen bei der initialisation und rebuild
	virtual G4VPhysicalVolume* Construct();//mit aufräumen

	//wird von Construct() bei initialisation und rebuild aufgerufen
	G4VPhysicalVolume* ConstructDetector();//ohne aufräumen

	//wird von Construct() bei initialisation aufgerufen.
	void InitializeMaterials();//einmalig ausgeführt.

	//wird vom Konstruktor aufgerufen
	void InitializeDimensions();//einmalig ausgeführt.

	//wird vom Konstruktor aufgerufen
	void InitializeRotations();//einmalig ausgeführt.

    	//Für die Hit-Generation über die eingebaute Funktionalität von G4VSensitiveDetector
    	//I think we rather distinguish our sensitive volume by hand in our
    	//L1kSteppingAction
	virtual void ConstructSDandField();

    	void UpdateGeometry();

    	//Messenger functions
    	//TODO

protected:
	G4VPhysicalVolume* worldPhys;
	G4VPhysicalVolume* cryostatPhys;
	G4VPhysicalVolume* fiberShroudInnerPhys;
	G4VPhysicalVolume* fiberShroudOuterPhys;
	G4VPhysicalVolume* lArPhys;
	G4VPhysicalVolume* wslrCopperPhys;
	G4VPhysicalVolume* wslrTetraTexPhys;
	G4VPhysicalVolume* wslrTPBPhys;

	//Materialien:
	G4Material* world_mat;
	G4Material* copper_mat;
	G4Material* lAr_mat;
	G4Material* tetraTex_mat;
	G4Material* steel_mat;
	G4Material* TPB_mat;

	//primäre Dimensionen (zwischenwerte werden in CostructDetector() angelegt und berechnet)
	//alle Längen sind NICHT halbiert
	G4double world_len, world_wid, world_height;//Längen Welt (in x, y, z)

	//For cryostat design
	G4double hneck;
   	G4double htopcylbot;
   	G4double hlittlecyl;
   	G4double slopetoplid;
   	G4double slopebotlid;

   	G4double rneck;
  	G4double rcyl;
   	G4double rlittlecyl;
   	G4double wallthickness;

	//for the inner fiber shroud
	G4double innerShroudInnerR;
	G4double innerShroudOuterR;
	G4double innerShroudHeight;

	//for the outer fiber shroud
	G4double outerShroudInnerR;
	G4double outerShroudOuterR;
	G4double outerShroudHeight;

	//for WSLR
	G4double wslrCopperInnerR;
	G4double wslrCopperOuterR;
	G4double wslrTetraTexInnerR;
	G4double wslrTetraTexOuterR;
	G4double wslrTPBThickness;
	G4double wslrHeight;



	void ConstructCryostat();
	void FillLAr();
	void BuildInnerShroud();
	void BuildOuterShroud();
	void BuildWSLRCopper();
	void BuildWSLRTetra();
	void BuildWSLRTPB();
	//void BuildWSLRChimney();
};

#endif
















