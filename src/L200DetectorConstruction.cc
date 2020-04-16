#include "L200DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Polycone.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4MaterialTable.hh"
#include "G4UnitsTable.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4VisAttributes.hh"
#include "G4Color.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RotationMatrix.hh"
#include "G4Region.hh"
#include "G4ios.hh"

#include "G4Transform3D.hh"

#include <cmath>

// = = = = = = = = = = = KONSTRUKTOR & DESTRUKTOR = = = = = = = = = = = = = =

L200DetectorConstruction::L200DetectorConstruction()
	: G4VUserDetectorConstruction()
{
	worldPhys 		= NULL;
	cryostatPhys 		= NULL;
	lArPhys			= NULL;
	fiberShroudInnerPhys 	= NULL;
	fiberShroudOuterPhys	= NULL;
	wslrCopperPhys		= NULL;
	wslrTetraTexPhys	= NULL;
	wslrTPBPhys		= NULL;

	world_mat 		= NULL;
	lAr_mat 		= NULL;
	steel_mat		= NULL;
	copper_mat		= NULL;
	tetraTex_mat		= NULL;
	TPB_mat			= NULL;


	InitializeDimensions();
	InitializeRotations();
}

L200DetectorConstruction::~L200DetectorConstruction()
{
}

//= = = = = = = = = = = = = = = INITIALISIERUNGEN = = = = = = = = = = = = = =

void L200DetectorConstruction::InitializeDimensions(){

	world_len = 10000*mm;
	world_wid = 10000*mm;
	world_height = 10000*mm;

	//cryostat values (GERDA values)
 	hneck = (1494. + 207.3) *mm+100.*mm;
   	htopcylbot = 5881.2*mm;
   	hlittlecyl = 221.0*mm;
   	slopetoplid = 0.24;
   	slopebotlid = 0.24;

   	rneck = 0.5*1000.0*mm;
   	rcyl = 2100.0*mm;
   	rlittlecyl =0.5*500.0*mm;

	wallthickness = 10*mm;

	//fiber shroud values
	innerShroudOuterR = 134.5*mm;
	innerShroudInnerR = 122.5*mm;
	innerShroudHeight = 1300*mm;

	outerShroudOuterR = 295*mm;
	outerShroudInnerR = 283*mm;
	outerShroudHeight = 1500*mm;

	//WSLR values
	wslrCopperOuterR = 310*mm; //TODO chcek dim
	wslrCopperInnerR = wslrCopperOuterR-0.03*mm;   //TODO check dim

	wslrTetraTexOuterR = wslrCopperInnerR;
	wslrTetraTexInnerR = wslrTetraTexOuterR- 0.01*mm; //TODO check dim

	wslrTPBThickness = 0.001*mm;
	wslrHeight = 3500*mm; //TODO check dim

}

void L200DetectorConstruction::InitializeRotations(){

}

void L200DetectorConstruction::InitializeMaterials(){
	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();
	G4double Z, M;

   	G4double density,fractionmass,abundance;
   	G4String name,symbol;
   	G4int ncomponents;
   	G4int natoms;
   	G4double temperature, pressure;
   	G4State  state;

	//Play chemestry
	G4Element* elH = new G4Element("H", "H", Z=1., M=1.01*g/mole);
	//G4Element* carbon = new G4Element("C", "C", Z=6., M=12.01*g/mole);
	//G4Element* oxygen = new G4Element("O", "O", Z=8., M=16.0*g/mole);
	//G4Element* silicon = new G4Element("Si", "Si", Z=14., M=28.09*g/mole);
	//G4Element* nitrogen = new G4Element("N", "N", Z=7., M=28.09*g/mole); //really this M?
	//G4Element* elLi = new G4Element(name="Lithium",symbol="Li",Z=3.,M=6.941*g/mole);
    	G4Element* argon = new G4Element("Argon","Ar",Z=18, M=39.948*g/mole);
	//G4Element* elCl = new G4Element(name="Chlorine",symbol="Cl",z=17.,a=35.45*g/mole);
	G4Element* elF = new G4Element(name="Fluorine","F",Z=9.,M=19.00*g/mole);
	G4Element* elC = new G4Element(name="Carbon","C",Z=6.,M=12.011*g/mole);

	//Some nice LAr
  	density=1.390*g/cm3;
  	ncomponents=1;
  	natoms=1;
  	state=kStateLiquid;
  	temperature=87.0*kelvin;
  	pressure=1.0*bar;
  	name="LiquidArgon";
   	lAr_mat = new G4Material(name,density,ncomponents,state,temperature,pressure);
        lAr_mat-> AddElement(argon, natoms = 1);

 	//Copper
	copper_mat = nist->FindOrBuildMaterial("G4_Cu");

	//Steel
	steel_mat =nist->FindOrBuildMaterial("G4_Stainless-Steel");

	//TetraTex = PTFE
	tetraTex_mat = new G4Material(name="tetraTex",
				      density=2.165*g/cm3,
				      ncomponents=2);
  	tetraTex_mat->AddElement(elF,fractionmass= 0.76);
  	tetraTex_mat->AddElement(elC,fractionmass= 0.24);

	//TPB
	TPB_mat = new G4Material(name="TPB",
				 density=1.08*g/cm3,
				 ncomponents=2);
	TPB_mat->AddElement(elC,22);
	TPB_mat->AddElement(elH,28);

	//World
	world_mat = nist->FindOrBuildMaterial("G4_Air");
}


// = = = = = = = = = = = = = = = CONSTRUCT = = = = = = = = = = = = = =

G4VPhysicalVolume* L200DetectorConstruction::Construct(){
	G4cout << "constructing..." << G4endl;
	//aufräumen, falls nötig:
	if (worldPhys != NULL) {
     		G4GeometryManager::GetInstance()->OpenGeometry();
     		G4PhysicalVolumeStore::GetInstance()->Clean();
     		G4LogicalVolumeStore::GetInstance()->Clean();
     		G4SolidStore::GetInstance()->Clean();
     		G4LogicalSkinSurface::CleanSurfaceTable();
     		G4LogicalBorderSurface::CleanSurfaceTable();
  	}

	if(world_mat == NULL){InitializeMaterials();}

	//(neu) aufbauen:
	return ConstructDetector();
}


G4VPhysicalVolume* L200DetectorConstruction::ConstructDetector(){

	//DEFINE VOLUMINA

	//World
  	G4Box* worldBox = new G4Box("World", world_len/2., world_wid/2.,world_height/2.);
  	G4LogicalVolume* worldLog = new G4LogicalVolume(worldBox, world_mat, "World");
  	worldPhys = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),worldLog, "World", 0, false, 0);

	//Place cryostat into the world
	G4cout << "Build cryostat" << G4endl;
	ConstructCryostat();
	G4cout << "Fill cryostat" << G4endl;
	FillLAr();
	G4cout << "Mount inner fiber shroud" << G4endl;
	BuildInnerShroud();
	G4cout << "Place outer fiber shroud" << G4endl;
	BuildOuterShroud();
	G4cout << "Lower WSLR" << G4endl;
	BuildWSLRCopper();
	BuildWSLRTetra();
	BuildWSLRTPB();

	// Instantiation of a set of visualization attributes with cyan colour
    	G4VisAttributes * cryostatVisAtt = new G4VisAttributes(G4Colour(0.,1.,1.)); //cyan
    	G4VisAttributes * detectorsVisAtt = new G4VisAttributes(G4Colour(1.,0.,1.)); //red
    	G4VisAttributes * worldVisAtt = new G4VisAttributes(G4Colour(0.,0.,0.)); //black

    	// Set the forced wireframe style
    	cryostatVisAtt->SetForceWireframe(true);
    	worldVisAtt->SetForceWireframe(true);


    	// Assignment of the visualization attributes to the logical volume
    	cryostatPhys->GetLogicalVolume()->SetVisAttributes(cryostatVisAtt);
    	worldPhys->GetLogicalVolume()->SetVisAttributes(worldVisAtt);

	return worldPhys;
}


//= = = = = = = = = = = = = = = SENSITIVE DETECTOR = = = = = = = = = = = = = =

void L200DetectorConstruction::ConstructSDandField(){

}


//The cryostat
void L200DetectorConstruction::ConstructCryostat()
{
	// Construct solid
    	G4double phistart(0.0 *deg);
    	G4double phitot(360.0 *deg);
    	G4int numzplanes(6);
    	G4double coord_z[6];
    	G4double coord_rInner[6] = {0,0,0,0,0,0};
    	G4double coord_rOuter[6] = {rneck,rneck,rcyl,rcyl,rlittlecyl,rlittlecyl};
    	//double coord_rOuter[6] = {rlittlecyl,rlittlecyl,rcyl,rcyl,rneck,rneck};

	for(int i=0; i<sizeof(coord_rInner)/sizeof(*coord_rInner); i++){
		coord_rInner[i] = coord_rOuter[i]-wallthickness;
	}

    	coord_z[0] = 0;
    	coord_z[1] = - hneck;
    	coord_z[2] = - hneck - slopetoplid*(rcyl-rneck);
    	coord_z[3] = - hneck - htopcylbot + slopebotlid*(rcyl-rlittlecyl);
    	coord_z[4] = - hneck - htopcylbot;
    	coord_z[5] = - hneck - htopcylbot - hlittlecyl;

    	G4double generalPolyconeZShift = -0.5*(coord_z[5] - hneck);

    	for(int i=0;i<numzplanes;i++) {
        	coord_z[i] += generalPolyconeZShift;
    	}

   	G4Polycone* thePolycone;
   	G4String theCryoPartName = "cryostat";

    	if(hneck==0) {
        	thePolycone = new G4Polycone(theCryoPartName,
        	                             phistart,phitot,
        	                             numzplanes-1,
               		                      &(coord_z[1]),&(coord_rInner[1]),&(coord_rOuter[1]));
    	}
   	 else {
        	thePolycone = new G4Polycone(theCryoPartName,
                	                     phistart,phitot,
                        	             numzplanes,
                                	     coord_z,coord_rInner,coord_rOuter);
    	}

    	//Logical Volume
    	G4LogicalVolume* cryostatLog = new G4LogicalVolume(thePolycone,steel_mat ,theCryoPartName);

    	//Physical volume
    	this->cryostatPhys = new G4PVPlacement(0,
						G4ThreeVector(0.,0.,0.),
						cryostatLog,"cryostat",
						this->worldPhys->GetLogicalVolume(),0,0);

}
//Liquid Argon
void L200DetectorConstruction::FillLAr(){
	// Construct solid
    	G4double phistart(0.0 *deg);
    	G4double phitot(360.0 *deg);
    	G4int numzplanes(6);
    	G4double coord_z[6];
    	G4double coord_rInner[6] = {0,0,0,0,0,0};
    	G4double coord_rOuter[6] = {rneck,rneck,rcyl,rcyl,rlittlecyl,rlittlecyl};
    	//double coord_rOuter[6] = {rlittlecyl,rlittlecyl,rcyl,rcyl,rneck,rneck};

	for(int i=0; i<sizeof(coord_rOuter)/sizeof(*coord_rOuter); i++){
		coord_rOuter[i] = coord_rOuter[i]-wallthickness;
	}

    	coord_z[0] = 0;
    	coord_z[1] = - hneck;
    	coord_z[2] = - hneck - slopetoplid*(rcyl-rneck);
    	coord_z[3] = - hneck - htopcylbot + slopebotlid*(rcyl-rlittlecyl);
    	coord_z[4] = - hneck - htopcylbot;
    	coord_z[5] = - hneck - htopcylbot - hlittlecyl;

    	G4double generalPolyconeZShift = -0.5*(coord_z[5] - hneck);

    	for(int i=0;i<numzplanes;i++) {
        	coord_z[i] += generalPolyconeZShift;
    	}

   	G4Polycone* thePolycone;
   	G4String theCryoPartName ="larVolume";

    	if(hneck==0) {
        	thePolycone = new G4Polycone(theCryoPartName,
        	                             phistart,phitot,
        	                             numzplanes-1,
               		                      &(coord_z[1]),&(coord_rInner[1]),&(coord_rOuter[1]));
    	}
   	 else {
        	thePolycone = new G4Polycone(theCryoPartName,
                	                     phistart,phitot,
                        	             numzplanes,
                                	     coord_z,coord_rInner,coord_rOuter);
    	}

    	//Logical Volume
    	G4LogicalVolume* lArLog = new G4LogicalVolume(thePolycone,lAr_mat ,theCryoPartName);

   	//Physical volume
    	this->lArPhys = new G4PVPlacement(0,
					  G4ThreeVector(0.,0.,0.),
					  lArLog,"larVolume",
					  this->cryostatPhys->GetLogicalVolume(),0,0);


}
//The inner fibershroud
void L200DetectorConstruction::BuildInnerShroud(){
	G4Tubs* isT = new G4Tubs("innerShroud",
				 innerShroudInnerR,
				 innerShroudOuterR,
				 innerShroudHeight/2.,
				 0,
				 2*M_PI);
	G4LogicalVolume* isLog = new G4LogicalVolume(isT,
						     this->lAr_mat,
						     "innerShroud");
	this->fiberShroudInnerPhys = new G4PVPlacement(0,
						G4ThreeVector(0.,0.,0.),
						isLog,
						"innerShroud",
						this->lArPhys->GetLogicalVolume(),0,0);
}

//The outer fibershroud
void L200DetectorConstruction::BuildOuterShroud(){
	G4Tubs* osT = new G4Tubs("outerShroud",
				 outerShroudInnerR,
				 outerShroudOuterR,
				 outerShroudHeight/2.,
				 0,
				 2*M_PI);
	G4LogicalVolume* osLog = new G4LogicalVolume(osT,
						     this->lAr_mat,
						     "outerShroud");
	this->fiberShroudOuterPhys = new G4PVPlacement(0,
						     G4ThreeVector(0.,0.,0.),
						     osLog,
						     "outerShroud",
						      this->lArPhys->GetLogicalVolume(),0,0);
}


//The WSLR Copper
void L200DetectorConstruction::BuildWSLRCopper(){
	G4Tubs* wslrcT = new G4Tubs("wslrCopper",
				 wslrCopperInnerR,
				 wslrCopperOuterR,
				 wslrHeight/2.,
				 0,
				 2*M_PI);
	G4LogicalVolume* wslrcLog = new G4LogicalVolume(wslrcT,
						     this->copper_mat,
						     "wslrCopper");
	this->wslrCopperPhys = new G4PVPlacement(0,
						G4ThreeVector(0.,0.,0.),
						wslrcLog,
						"wslrCopper",
						this->lArPhys->GetLogicalVolume(),0,0);

}

//The WSLR TetraTex
void L200DetectorConstruction::BuildWSLRTetra(){
	G4Tubs* wslrtT = new G4Tubs("wslrTetraTex",
				 wslrTetraTexInnerR,
				 wslrTetraTexOuterR,
				 wslrHeight/2.,
				 0,
				 2*M_PI);
	G4LogicalVolume* wslrtLog = new G4LogicalVolume(wslrtT,
						     this->tetraTex_mat,
						     "wslrTetra");
	this->wslrTetraTexPhys = new G4PVPlacement(0,
						G4ThreeVector(0.,0.,0.),
						wslrtLog,
						"wslrTetra",
						this->lArPhys->GetLogicalVolume(),0,0);

}

//The WSLR TPB
void L200DetectorConstruction::BuildWSLRTPB(){
	G4Tubs* wslrTPBT = new G4Tubs("wslrTPB",
				 wslrTetraTexInnerR-wslrTPBThickness,
				 wslrTetraTexInnerR,
				 wslrHeight/2.,
				 0,
				 2*M_PI);
	G4LogicalVolume* wslrTPBLog = new G4LogicalVolume(wslrTPBT,
						     this->TPB_mat,
						     "wslrTPB");
	this->wslrTPBPhys = new G4PVPlacement(0,
					   G4ThreeVector(0.,0.,0.),
					   wslrTPBLog,
					   "wslrTPB",
					   this->lArPhys->GetLogicalVolume(),0,0);

}

void L200DetectorConstruction::UpdateGeometry()
{
	G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}







