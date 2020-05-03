#include "Randomize.hh"

#include "G4ios.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Alpha.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleGun.hh"
#include "G4Run.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"

#include "L200ParticleGenerator.hh"
#include "L200ParticleGeneratorMessenger.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

using namespace CLHEP;

const G4double L200ParticleGenerator::LambdaE = twopi *1.973269602e-16 * m * GeV;

L200ParticleGenerator::L200ParticleGenerator()
	: scanAngle(2*M_PI/28), flatVoxelIndex(0), verbosity(0), abortVoxel(false), abortOnNonlar(true)
{
	fMessenger = new L200ParticleGeneratorMessenger(this);
	fParticleGun = new G4ParticleGun(1);
	fLArWL = 128*nm;
	fCurrentEnergy = LambdaE/fLArWL;
}


L200ParticleGenerator::~L200ParticleGenerator()
{
  delete fMessenger;
  delete fParticleGun;
}


int L200ParticleGenerator::nextVoxel(){

	abortVoxel = false;

	//#### Part I: make voxel pattern over 1st quadrant ###
	G4double xMax = fRadiusMax;
	G4double xMin = 0.;//-fRadiusMax;
  	G4int xBins = (xMax-xMin)/fBinWidth;///fRadiusMax*2/fBinWidth;
  	G4double yMax  = fRadiusMax;
	G4double yMin = 0.;//-fRadiusMax;
  	G4int yBins = (yMax-yMin)/fBinWidth;//fRadiusMax*2/fBinWidth;
	G4int zBins = 1;
	G4double zMin = fZ;
	G4double zMax = fZ;

	switch(fDim){
	case 1:
		yMax = 0;
		yMin = 0;
		yBins = 1;
		xMin = 0;
		xBins = fRadiusMax/fBinWidth;
		break;
	case 2:
		break;
	case 3:
		zMin =  0;
		zBins = (zMax-zMin)/fBinWidth;
		break;
	default:
		break;
	}

  if(flatVoxelIndex == 0 && verbosity >= 1){
    G4cout<<"N bins XY "<<xBins*yBins<<" x/y Max "<<xMax/cm<<"(cm), x/y Min "<<fRadiusMin/cm<<"(cm), Z "<<fZ<<"(cm), binWidth "<<fBinWidth/cm<<"(cm)..."<< fNParticles<<" particles per voxel, with "<<fNParticles*xBins*yBins<<" photons generated"<<G4endl;
  }

	//### Part II: define dimensions of current voxel; skipping if exceeds angle ###
	while(true){	//loop should not affect 1D case (break always after 1st call)
		if(flatVoxelIndex == xBins*yBins*zBins) return 0;		//escape condition

		//decide indices
		int iBinZ = flatVoxelIndex / (xBins*yBins);		//int division @ wörk
		int iBinY = (flatVoxelIndex-iBinZ*xBins*yBins)/xBins;
		int iBinX = flatVoxelIndex-iBinZ*xBins*yBins - iBinY*xBins;

		//decide positions
		  G4double x = xMin + iBinX*fBinWidth;
		  G4double y = yMin + iBinY*fBinWidth;
		  G4double z = zMin + iBinZ*fBinWidth;
		  //G4cout<<"Generating point randomly ("<<x<<","<<y<<","<<fZ<<")"<<G4endl;

		currentVoxel.xPos = x + fCenterVector.getX();
		currentVoxel.yPos = y + fCenterVector.getY();
		currentVoxel.zPos = z + fCenterVector.getZ();

		currentVoxel.xWid = fBinWidth;
		currentVoxel.yWid = fBinWidth;
		currentVoxel.zWid = fBinWidth;

		//escape skip loop in case bottom right point of voxel within angle
		// and bottom left point is still within radius
		if(currentVoxel.yPos <= tan(scanAngle)*(currentVoxel.xPos+currentVoxel.xWid) &&
			currentVoxel.xPos*currentVoxel.xPos+currentVoxel.yPos*currentVoxel.yPos <= fRadiusMax*fRadiusMax) break;

		if(verbosity >= 3) G4cout << "Skipping voxel "<<currentVoxel<<" for symmetry reasons" << G4endl;

		flatVoxelIndex++;	//increment if volume skipped
	}

	if(verbosity >= 3) G4cout << "Will use voxel "<<currentVoxel<<" in the following." << G4endl;

	//### Part III: increment index for next call & report particle count ###
	flatVoxelIndex++;	//increment after --> 1st voxel is 0
	return fNParticles;
}











void L200ParticleGenerator::DirectionDecider()
{
  G4double phi = 2*pi*G4UniformRand();
  G4double costheta = 2*G4UniformRand() -1;
  G4double theta = acos(costheta);

  G4double px = cos(phi)*sin(theta);
  G4double py = sin(phi)*sin(theta);
  G4double pz = cos(theta);

  fDirection.setX(px);fDirection.setY(py);fDirection.setZ(pz);
}




void L200ParticleGenerator::PositionDecider(G4double xPos,G4double yPos,G4double zPos,G4double binWidth)
{

  larFailed = false;
  G4ThreeVector rpos(1,1,1);
  G4bool isIn = false;
  int errorCounter = 0;
  while(!isIn){
	//Random position in voxel
  	G4double x = xPos,y=yPos,z=zPos;
  	G4double rand = G4UniformRand();
  	x += binWidth*rand;
  	rand = G4UniformRand();
  	y += binWidth*rand;
  	rand = G4UniformRand();
  	z += binWidth*rand;
  	rpos.setX(x);rpos.setY(y);rpos.setZ(z);

  	//Is it in the Argon?
  	//Note that in the Stepping Action, a step or two are taken before intial position is stored
  	isIn = IsInArgon(rpos);
  	errorCounter++;
  	if((!isIn) && (errorCounter > 1e2 || abortOnNonlar)){
    	   rpos.setX(1000000);rpos.setY(1000000);rpos.setZ(1000000);
	   // if error counter exceeds, we do not want to produce a photon hence the large values of co-ordinates
			larFailed = true;
    	   if(verbosity >= 0) G4cout<< "LAr error has exceeded 100" << G4endl;
    	   break;
  	}
  }
  if(verbosity >= 4) G4cout<<"Generator vertex "<<rpos<<G4endl;
  fCurrentPosition = rpos;

}
G4bool L200ParticleGenerator::IsInArgon(G4ThreeVector rpos)
{
  bool isit = false;
  //This is how Geant4 suggests you should randomly generate a point in a volume
  G4ThreeVector myPoint = rpos;
  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4VPhysicalVolume* myVolume = theNavigator->LocateGlobalPointAndSetup(myPoint);
  if(myVolume->GetName() == "larVolume") isit = true;
  //G4cout<< " The current material is " << myVolume->GetName() << " okay " <<G4endl;
  return isit;
}


void L200ParticleGenerator::GeneratePrimaryVertex(G4Event *event)
{
	if(abortVoxel) return;		//dont mess around any more with a aborted voxel.

    fParticleGun->SetParticlePolarization(G4ThreeVector(2*G4UniformRand()-1,2*G4UniformRand()-1,2*G4UniformRand()-1 ) );

    //what is the particle
    fParticleGun->SetParticleDefinition(G4OpticalPhoton::OpticalPhotonDefinition());
    //determine particle momentum direction
    DirectionDecider();

    //determine particle position
    PositionDecider(currentVoxel.xPos,
		    currentVoxel.yPos,
		    currentVoxel.zPos,
		    fBinWidth);
    //if(fCurrentPosition == G4ThreeVector(1000000,1000000,1000000))return;// break;
	if(larFailed){		//fail bit arrived from position decider
		if(abortOnNonlar){
			abortVoxel = true;
			if(verbosity >= 1) G4cout << "Aborting voxel "<<currentVoxel<<G4endl;
		}
		return;		//go out of primary production immediately
	}

    //particle direction, position, and energy sent to ParticleGun
    fParticleGun->SetParticlePosition(fCurrentPosition);
    fParticleGun->SetParticleMomentumDirection(fDirection);
    fParticleGun->SetParticleEnergy(fCurrentEnergy);
    fParticleGun->SetNumberOfParticles(1);

    //vertex generated by ParticleGun
    fParticleGun->GeneratePrimaryVertex(event);
}













