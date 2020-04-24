#include "L200DetectorConstruction.hh"
#include "L200DetectorMessenger.hh"

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
#include "G4SurfaceProperty.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RotationMatrix.hh"
#include "G4Region.hh"
#include "G4ios.hh"
#include "G4AssemblyVolume.hh"


#include "G4Transform3D.hh"

#include <cmath>
using namespace CLHEP;

// = = = = = = = = = = = KONSTRUKTOR & DESTRUKTOR = = = = = = = = = = = = = =

L200DetectorConstruction::L200DetectorConstruction()
	: G4VUserDetectorConstruction(), geAssembly(NULL)
{
	det_briefTaube = new L200DetectorMessenger(this);

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
	enrGe_mat = NULL;

	InitializeDimensions();
	InitializeRotations();
}

L200DetectorConstruction::~L200DetectorConstruction()
{
	delete det_briefTaube;
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
	wlsrRadius = 700*mm;
	wlsrCuThickness = 0.03*mm;

	wlsrTetraTexThickness = 0.01*mm;

	wslrTPBThickness = 0.001*mm;
	wslrHeight = 3500*mm;

	geDiscHeight = 100*mm;
	geDiscRad = 50*mm;
	geDiscGap = 50*mm;
	geArrayRad = 200*mm;			//TODO finalize at some point
	geDetectorsInString = 7;		//should be 7*10cm+6*5cm = 1m total
	geStringCount = 14;				//defines the symmetry angle: 360/(14*2) ~ 12.8°

	//LAr optical
	lArAbsVUV = 20*cm;
	lArAbsVis = 1000*m;
	lArWL = 128*nm;
	tpbWL = 450*nm;
	deltaE = 0.1*electronvolt;	//half width of emission / absorption peaks
	lambdaE = twopi*1.97326902e-16 * m * GeV;		//hc/e
	lArRay = false;
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

	G4Isotope* Ge70 = new G4Isotope(name="Ge70",  32, 70, 69.92*g/mole);
  	G4Isotope* Ge72 = new G4Isotope(name="Ge72",  32, 72, 71.92*g/mole);
  	G4Isotope* Ge73 = new G4Isotope(name="Ge73",  32, 73, 73.0*g/mole);
 	G4Isotope* Ge74 = new G4Isotope(name="Ge74",  32, 74, 74.0*g/mole);
  	G4Isotope* Ge76 = new G4Isotope(name="Ge76",  32, 76, 76.0*g/mole);

	M = 75.71*g/mole;
  	G4int nIsotopes=5;
	G4Element* elGeEnr = new G4Element(name="enrichedGermanium", symbol="GeEnr",nIsotopes);
  	elGeEnr->AddIsotope(Ge70,abundance= 0.0*perCent);
  	elGeEnr->AddIsotope(Ge72,abundance= 0.1*perCent);
  	elGeEnr->AddIsotope(Ge73,abundance= 0.2*perCent);
  	elGeEnr->AddIsotope(Ge74,abundance= 13.1*perCent);
  	elGeEnr->AddIsotope(Ge76,abundance= 86.6*perCent);

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
	//optical LAr

	const G4int NUM = 2;
	const G4int NUM_lines = 6;		//for emission/absorption --> no interpolation, but deltas


	G4double photonEnergy[NUM] = {lambdaE/(lArWL),lambdaE/(tpbWL)};		//defines energy from large to small
	G4cout << "Photonenergy";
	G4cout << "Photon energy of " << lArWL/nm << " nm is: " <<photonEnergy[0]/eV << " eV" << G4endl;
	G4double lArAbsorption[NUM] = {lArAbsVUV,lArAbsVis};
	G4double lArRefIndex[NUM]={LArRefIndex(lArWL),LArRefIndex(tpbWL)};
	G4double lArRayLength[NUM]={LArRayLength(lArWL,temperature),LArRayLength(tpbWL,temperature)};

	G4MaterialPropertiesTable* mptLAr = new G4MaterialPropertiesTable();
	mptLAr->AddProperty("RINDEX",photonEnergy,lArRefIndex,NUM);
	if(lArRay)
		mptLAr->AddProperty("RAYLEIGH",photonEnergy,lArRayLength,NUM);

	mptLAr->AddProperty("ABSLENGTH",photonEnergy,lArAbsorption,NUM);

	lAr_mat->SetMaterialPropertiesTable(mptLAr);

 	//Copper
	copper_mat = nist->FindOrBuildMaterial("G4_Cu");

	//Steel
	steel_mat =nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");

	//TetraTex = PTFE
	tetraTex_mat = new G4Material(name="tetraTex",
				      density=2.165*g/cm3,
				      ncomponents=2);
  	tetraTex_mat->AddElement(elF,fractionmass= 0.76);
  	tetraTex_mat->AddElement(elC,fractionmass= 0.24);

	//TPB
	G4double photonEnergy_lines[NUM_lines] = {
		//9.8*eV, 9.7*eV, 9.6*eV, 2.85*eV, 2.75*eV, 2.65*eV};
		//2.65*eV, 2.75*eV, 2.85*eV, 8.6*eV, 9.7*eV, 10.8*eV};
	//remember: energy goes from large to small (here)
			lambdaE/(lArWL) + deltaE, lambdaE/(lArWL), lambdaE/(lArWL) - deltaE,
			lambdaE/(tpbWL) + deltaE, lambdaE/(tpbWL), lambdaE/(tpbWL) - deltaE};
	TPB_mat = new G4Material(name="TPB",
				 density=1.08*g/cm3,
				 ncomponents=2);
	TPB_mat->AddElement(elC,22);
	TPB_mat->AddElement(elH,28);
	//optical TPB
	G4double tpbQuantumEff = 1.2;
	G4double tpbTimeConst = 0.01*ns;
	G4double tpbRefIndex[NUM_lines] = {1.635,1.635,1.635,1.635,1.635,1.635};
	G4double tpbEmission[NUM_lines] = {0., 0., 0., 0., 1., 0.};//{lambdaE/tpbWL,lambdaE/tpbWL};
	G4double tpbAbsorption[NUM_lines] = {1*nm, 1*nm, 1*nm, 1000*m, 1000*m, 1000*m};
//{1000*m, 1000*m, 10*nm, 1*nm, 1*nm, 1*nm};
	G4MaterialPropertiesTable* tbpMPT = new G4MaterialPropertiesTable();
	tbpMPT->AddProperty("RINDEX",photonEnergy_lines,tpbRefIndex,NUM_lines);
	tbpMPT->AddProperty("WLSABSLENGTH",photonEnergy_lines, tpbAbsorption,NUM_lines);
	tbpMPT->AddProperty("WLSCOMPONENT",photonEnergy_lines, tpbEmission, NUM_lines);
	tbpMPT->AddConstProperty("WLSTIMECONSTANT", tpbTimeConst);
	tbpMPT->AddConstProperty("WLSMEANNUMBERPHOTONS", tpbQuantumEff);
	TPB_mat->SetMaterialPropertiesTable(tbpMPT);

	
	// enriched germanium
  	density = 5.56*g/cm3;//cryogenic temperature (at room temperature: 5.54*g/cm3)
  	enrGe_mat = new G4Material(name="EnrichedGe", density, 1);
  	enrGe_mat->AddElement(elGeEnr,natoms=1);

	//World
	world_mat = nist->FindOrBuildMaterial("G4_AIR");

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
	G4cout << "Now come the detectors!"<<G4endl;
	BuildGeDetectors();
	G4cout << "Build optics" << G4endl;
	BuildOptics();

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

	G4cout << "Doing a sanity check. Keep fingers crossed!"<<G4endl;
	sanityCheck();

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

	/*for(int i=0; i<sizeof(coord_rInner)/sizeof(*coord_rInner); i++){
		coord_rInner[i] = coord_rOuter[i]-wallthickness;
	}*/

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

	//optical stuff



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
	wslrCopperOuterR = wlsrRadius;
	wslrCopperInnerR = wslrCopperOuterR-wlsrCuThickness;
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
	wslrTetraTexOuterR = wslrCopperInnerR;
	wslrTetraTexInnerR = wslrTetraTexOuterR- wlsrTetraTexThickness;
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

//detectörs
//without any optical surface added to the detectors, I see always LambertianReflection of 128 nm photons.
//I guess it's due to the LocgicalSkinSurface attached to the larVolume
//lets see if it gets overwritten by a LogicalBorderSurface...
void L200DetectorConstruction::BuildGeDetectors(){
	G4Tubs* geDisc_tub = new G4Tubs("Ge_detector",0,geDiscRad,geDiscHeight/2.,0,2*M_PI);
	/*G4LogicalVolume* */geDisc_log = new G4LogicalVolume(geDisc_tub,enrGe_mat,"Ge_detector");
	
	G4double stringFullHeight = geDetectorsInString*geDiscHeight+(geDetectorsInString-1)*geDiscGap;

	//using again the assembly to create a nice pattern
	delete geAssembly;	//delete old one to prevent leak; deleting NULL is unproblematic
	/*G4AssemblyVolume* */geAssembly = new G4AssemblyVolume();
	//location of disc within assembly:
	G4RotationMatrix Ra(0.,0.,0.); G4ThreeVector radialVect; 
	//location of assembly in LAr:
	G4RotationMatrix Rm = G4RotationMatrix(0.,0.,0.); G4ThreeVector Tm = G4ThreeVector(0.,0.,0.);

	for(int i = 0; i < geStringCount; i++){		//outer: different strings
		radialVect.setRThetaPhi(geArrayRad,0.5*M_PI,2.*i*M_PI/geStringCount);
		for(int j = 0; j < geDetectorsInString; j++){
			G4ThreeVector verticalVect(0.,0.,-stringFullHeight/2.+(j+0.5)*geDiscHeight+j*geDiscGap);
			G4ThreeVector sum = radialVect + verticalVect;
			geAssembly->AddPlacedVolume(geDisc_log, sum, &Ra);
		}
	}

	geAssembly->MakeImprint(lArPhys->GetLogicalVolume(), Tm, &Rm);
	
}


void L200DetectorConstruction::UpdateGeometry()
{
	G4cout << "Geometry updated" << G4endl;
	G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}


void L200DetectorConstruction::BuildOptics()
{
	//TPB <-> LAr
    	G4OpticalSurface* osIn = new G4OpticalSurface("LArToTPB",unified,ground,dielectric_dielectric,.5);
	const G4int NUM = 2;
	G4double photonEnergy[NUM] = {lambdaE/(tpbWL), lambdaE/(lArWL)};
	G4double tpbRefIndex[NUM] = {1.635,1.635};
	G4MaterialPropertiesTable *mptInTPB = new G4MaterialPropertiesTable();
	mptInTPB->AddProperty("RINDEX",photonEnergy,tpbRefIndex,NUM);
	osIn->SetMaterialPropertiesTable(mptInTPB);

	new G4LogicalBorderSurface("LAr_TO_WLSRTPB",lArPhys,wslrTPBPhys,osIn);
	new G4LogicalBorderSurface("WLSRTPB_TO_LAr",wslrTPBPhys,lArPhys,osIn);

	//Tetratex
	G4OpticalSurface* sf = new G4OpticalSurface("TetraTex_Surface",unified,groundfrontpainted, dielectric_dielectric);
	G4double tetraReflectivity[NUM] = {0.95,0.95};
	G4MaterialPropertiesTable *mptIntetra = new G4MaterialPropertiesTable();
	mptIntetra->AddProperty("REFLECTIVITY",photonEnergy,tetraReflectivity,NUM);
	sf->SetMaterialPropertiesTable(mptIntetra);
	new G4LogicalSkinSurface("TetraTex_Surface",wslrTetraTexPhys->GetLogicalVolume(),sf);


	//optical stuff for cu
	G4OpticalSurface* sfcu = new G4OpticalSurface("Cu_surface",unified,ground,dielectric_metal,0.5);
	G4double cuReflectivity[NUM] = {0.4448,0.15};
	G4double cuEfficiency[NUM] = {0.,0.};
	G4MaterialPropertiesTable *mptInCu = new G4MaterialPropertiesTable();
	mptInCu->AddProperty("REFLECTIVITY",photonEnergy,cuReflectivity,NUM);
	mptInCu->AddProperty("EFFICIENCY",photonEnergy,cuEfficiency,NUM);
	sfcu->SetMaterialPropertiesTable(mptInCu);
	new G4LogicalSkinSurface("Cu_Surface",wslrCopperPhys->GetLogicalVolume(),sfcu);

	//optical stuff inner shroud
	//If the photon goes into the fiber, absorpe it
	G4OpticalSurface* sfIn = new G4OpticalSurface("LAr_TO_InnerFiber",unified,ground,dielectric_dielectric);
	G4double fiberReflectivity[NUM] = {0.,0.};
	G4double fiberEfficiency[NUM] = {1.,1.};
	G4MaterialPropertiesTable *mptIn = new G4MaterialPropertiesTable();
	mptIn->AddProperty("REFLECTIVITY",photonEnergy,fiberReflectivity,NUM);
	mptIn->AddProperty("EFFICIENCY",photonEnergy,fiberEfficiency,NUM);
	sfIn->SetMaterialPropertiesTable(mptIn);
	new G4LogicalBorderSurface("LAr_TO_InnerFiber",lArPhys,fiberShroudInnerPhys,sfIn);

	G4OpticalSurface* sfOut = new G4OpticalSurface("InnerFiber_TO_LAr",unified,ground,dielectric_dielectric);
	G4double fiberReflectivity2[NUM] = {0.,0.};
	G4double fiberEfficiency2[NUM] = {1.,1.};
	G4MaterialPropertiesTable *mptOut = new G4MaterialPropertiesTable();
	mptOut->AddProperty("REFLECTIVITY",photonEnergy,fiberReflectivity2,NUM);
	mptOut->AddProperty("EFFICIENCY",photonEnergy,fiberEfficiency2,NUM);
	sfIn->SetMaterialPropertiesTable(mptOut);
	new G4LogicalBorderSurface("InnerFiber_TO_LAr",fiberShroudInnerPhys,lArPhys,sfOut);

	//optical stuff outer shroud

	//If the photon goes into the fiber, absorpe it
	G4OpticalSurface* sfInOuter = new G4OpticalSurface("LAr_TO_OuterFiber",unified,ground,dielectric_dielectric);
	G4double outerfiberReflectivity[NUM] = {0.,0.};
	G4double outerfiberEfficiency[NUM] = {1.,1.};
	G4MaterialPropertiesTable *mptInOuter = new G4MaterialPropertiesTable();
	mptInOuter->AddProperty("REFLECTIVITY",photonEnergy,outerfiberReflectivity,NUM);
	mptInOuter->AddProperty("EFFICIENCY",photonEnergy,outerfiberEfficiency,NUM);
	sfInOuter->SetMaterialPropertiesTable(mptInOuter);
	new G4LogicalBorderSurface("LAr_TO_OuterFiber",lArPhys,fiberShroudOuterPhys,sfInOuter);

	G4OpticalSurface* sfOutOuter = new G4OpticalSurface("OuterFiber_TO_LAr",unified,ground,dielectric_dielectric);
	G4double outerfiberReflectivity2[NUM] = {0.,0.};
	G4double outerfiberEfficiency2[NUM] = {1.,1.};
	G4MaterialPropertiesTable *mptOutOuter = new G4MaterialPropertiesTable();
	mptOutOuter->AddProperty("REFLECTIVITY",photonEnergy,outerfiberReflectivity2,NUM);
	mptOutOuter->AddProperty("EFFICIENCY",photonEnergy,outerfiberEfficiency2,NUM);
	sfOutOuter->SetMaterialPropertiesTable(mptOutOuter);
	new G4LogicalBorderSurface("OuterFiber_TO_LAr",fiberShroudOuterPhys,lArPhys,sfOutOuter);

	//ge optical surface (only inward; light coming out of Ge should be rare...)
	//taken from MaGe: /legendgeometry/src/LGND_200_OpticalSurfaces.cc | 275 ("LArToGe")	
	// and materials/src/MGLGNDOpticalMaterialProperties.cc
	//reflectivity from generators/data/Reflectivity_Ge.dat (copied to data folder)	
	{		//using nice empty scope to prevent local name clash
	G4OpticalSurface* sfGe = new G4OpticalSurface("LArToGe",unified,groundfrontpainted,dielectric_metal,0.5);	//last val: smoothness
	G4double reflectivity[NUM] = {0.3563,0.6500};//keep in mind: optical (450nm), VUV
	G4double absorption[NUM] = {1*nm,1*nm};	//dummy; from MaGe
	G4double rIndex[NUM] = {1.25, 1.25};	//like in MaGe; commented out in MPT, however
	G4double efficiency[NUM] = {0., 0.};	//unlike MaGe, we don't think that Ge detects light very well
	G4MaterialPropertiesTable *mpt = new G4MaterialPropertiesTable();
	mpt->AddProperty("REFLECTIVITY",photonEnergy,reflectivity,NUM);
	mpt->AddProperty("EFFICIENCY",photonEnergy,efficiency,NUM);
	mpt->AddProperty("ABSLENGTH", photonEnergy, absorption, NUM);
  //mpt->AddProperty("RINDEX", photonEnergy, rIndex, NUM);
	sfGe->SetMaterialPropertiesTable(mpt);
	new G4LogicalSkinSurface("Ge_Surface",geDisc_log,sfGe);
	}
	

	//Give LAr volume a refractive skin
	G4double lArRefIndex[NUM]={LArRefIndex(tpbWL),LArRefIndex(lArWL)};
	G4MaterialPropertiesTable *mptLAr = new G4MaterialPropertiesTable();
	mptLAr->AddProperty("RINDEX",photonEnergy,lArRefIndex,NUM);

	G4OpticalSurface* sfLAr = new G4OpticalSurface("LAr_Surface",unified,groundfrontpainted, dielectric_dielectric);
	sfLAr->SetMaterialPropertiesTable(mptLAr);
	new G4LogicalSkinSurface("LAr_Surface",lArPhys->GetLogicalVolume(),sfLAr);

	



}

G4double L200DetectorConstruction::LArRefIndex(G4double lambda)
{
	G4cout << "LAr refindex for ";
	G4cout << lambda;
	G4cout << " is: ";
	G4double ret = sqrt(LArEpsilon(lambda));
  	G4cout << ret << G4endl;
	return ret; // square root of dielectric constant
}

// Calculates the dielectric constant of LAr from the Bideau-Sellmeier formula.
// See : A. Bideau-Mehu et al., "Measurement of refractive indices of Ne, Ar,
// Kr and Xe ...", J. Quant. Spectrosc. Radiat. Transfer, Vol. 25 (1981), 395

G4double L200DetectorConstruction::LArEpsilon(G4double lambda)
{
  G4double epsilon;
  if (lambda < 110*nanometer) return 1.0e4; // lambda MUST be > 110.0 nm
  epsilon = lambda / micrometer; // switch to micrometers
  epsilon = 1.0 / (epsilon * epsilon); // 1 / (lambda)^2
  epsilon = 1.2055e-2 * ( 0.2075 / (91.012 - epsilon) +
                          0.0415 / (87.892 - epsilon) +
                          4.3330 / (214.02 - epsilon) );
  epsilon *= (8./12.); // Bideau-Sellmeier -> Clausius-Mossotti
  G4double LArRho = 1.396*g/cm3;
  G4double GArRho = 1.66e-03*g/cm3;
  epsilon *= (LArRho / GArRho); // density correction (Ar gas -> LAr liquid)
  if (epsilon < 0.0 || epsilon > 0.999999) return 4.0e6;
  epsilon = (1.0 + 2.0 * epsilon) / (1.0 - epsilon); // solve Clausius-Mossotti
  return epsilon;
}

//-------------------------------------------------------------------------><
// Calculates the Rayleigh scattering length using equations given in
// G. M. Seidel at al., "Rayleigh scattering in rare-gas liquids",
// arXiv:hep-ex/0111054 v2 22 Apr 2002

G4double L200DetectorConstruction::LArRayLength(G4double lambda, G4double temp)
{
  G4double dyne = 1.0e-5*newton;
  static const G4double LArKT = 2.18e-10 * cm2/dyne; // LAr isothermal compressibility
  static const G4double k = 1.380658e-23 * joule/kelvin; // the Boltzmann constant
  G4double h;
  h = LArEpsilon(lambda);
  if (h < 1.00000001) h = 1.00000001; // just a precaution
  h = (h - 1.0) * (h + 2.0); // the "dielectric constant" dependance
  h *= h; // take the square
  h *= LArKT * temp * k; // compressibility * temp * Boltzmann constant
  h /= lambda * lambda * lambda * lambda; // (lambda)^4
  h *= 9.18704494231105429; // (2 * Pi / 3)^3
  if ( h < (1.0 / (10.0 * km)) ) h = 1.0 / (10.0 * km); // just a precaution
  if ( h > (1.0 / (0.1 * nanometer)) ) h = 1.0 / (0.1 * nanometer); // just a precaution
  return ( 1.0 / h );
}



//metod with licence to kill the run
void L200DetectorConstruction::sanityCheck(){
	if(hneck+htopcylbot+hlittlecyl>=world_height || rcyl >= world_len || rcyl >= world_wid ){
		G4Exception("L200DetectorConstruction::sanityCheck","cryoCrashesWorld",FatalException,"volume dimension conflict between world & cryostat (i.e. cryostat is too fat)");
	}
	if(geDiscHeight*geDetectorsInString + geDiscGap*(geDetectorsInString-1) > htopcylbot ){
		G4Exception("L200DetectorConstruction::sanityCheck","stringsTooLong",FatalException,"ge detector strings exceed cryo height");
	}
	if(geArrayRad - geDiscRad <= innerShroudOuterR){
		G4Exception("L200DetectorConstruction::sanityCheck", "stringTouchInnerShroud",FatalException,"ge detector strings touch inner fiber shroud");	
	}
	if(geArrayRad + geDiscRad >= outerShroudInnerR){
		G4Exception("L200DetectorConstruction::sanityCheck", "stringTouchOuterShroud",FatalException,"ge detector strings touch outer fiber shroud");	
	}
	//TODO further checks
	

	G4cout << "Sanity check finished successful." << G4endl;
}


























