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
	: flatVoxelIndex(0)
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
			

	G4double xMax = fRadiusMax, xMin = -fRadiusMax;
  G4int xBins = fRadiusMax*2/fBinWidth;
  G4double yMax  = fRadiusMax, yMin = -fRadiusMax;
  G4int yBins = fRadiusMax*2/fBinWidth;
  if(fis1D){
    yMax = 0;
    yMin = 0;
    yBins = 1;
    xMin = 0;
    xBins = fRadiusMax/fBinWidth;
  }
  if(flatVoxelIndex == 0){
    G4cout<<"N bins XY "<<xBins*yBins<<" x/y Max "<<xMax/cm<<"(cm), x/y Min "<<fRadiusMin/cm<<"(cm), Z "<<fZ<<"(cm), binWidth "<<fBinWidth/cm<<"(cm)..."<< fNParticles<<" particles per voxel, with "<<fNParticles*xBins*yBins<<" photons generated"<<G4endl;
  }
	if(flatVoxelIndex == xBins*yBins) return 0;		//escape condition

	//decide indices
	int iBinY = flatVoxelIndex / xBins;		//int division @ wörk
	int iBinX = flatVoxelIndex % xBins;


	//decide positions
      G4double x = xMin + iBinX*fBinWidth;
      G4double y = yMin + iBinY*fBinWidth;
      G4cout<<"Generating point randomly ("<<x<<","<<y<<","<<fZ<<")"<<G4endl;

	currentVoxel.xPos = x;
	currentVoxel.yPos = y;
	currentVoxel.zPos = fZ;				//TODO here if z free

	currentVoxel.xWid = fBinWidth;
	currentVoxel.yWid = fBinWidth;			//could y also make 0. width
	currentVoxel.zWid = 0.;
	



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
  	if(errorCounter > 1e2){
    	   rpos.setX(1000000);rpos.setY(1000000);rpos.setZ(1000000);
	   // if error counter exceeds, we do not want to produce a photon hence the large values of co-ordinates
    	   fCurrentPosition = rpos;
    	   G4cout<< "LAr error has exceeded 100" << G4endl;
    	   break;
  	}
  }
  G4cout<<"Generator vertex "<<rpos<<G4endl;
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
  G4cout<< " The current material is " << myVolume->GetName() << " okay " <<G4endl;
  return isit;
}


void L200ParticleGenerator::GeneratePrimaryVertex(G4Event *event)
{
  //User inputs specific height, radius, bin size, and number of points per bin in macro
  //generator creates N vertices at random postion in each voxel (defined by fBinWidth)
  //A cube of size fRadius at height fZ is scanned voxel by voxel
  //This method gaurentees a uniform number of points in each voxel
/*
  G4double xMax = fRadiusMax, xMin = -fRadiusMax;
  G4double xBins = fRadiusMax*2/fBinWidth;
  G4double yMax  = fRadiusMax, yMin = -fRadiusMax;
  G4double yBins = fRadiusMax*2/fBinWidth;
  if(fis1D){
    yMax = 0;
    yMin = 0;
    yBins = 1;
    xMin = 0;
    xBins = fRadiusMax/fBinWidth;
  }
  if(event->GetEventID() == 0){
    G4cout<<"N bins XY "<<xBins*yBins<<" x/y Max "<<xMax/cm<<"(cm), x/y Min "<<fRadiusMin/cm<<"(cm), Z "<<fZ<<"(cm), binWidth "<<fBinWidth/cm<<"(cm)..."<< fNParticles<<" particles per voxel, with "<<fNParticles*xBins*yBins<<" photons generated"<<G4endl;
  }
  for(int k = 0; k < yBins;k++){
    for(int j = 0; j < xBins;j++){
      G4double x;
      x = xMin + j*fBinWidth;
      G4double y;
      y = yMin + k*fBinWidth;
      G4cout<<"Generating point randomly ("<<x<<","<<y<<","<<fZ<<")"<<G4endl;*/

      //for(int i = 0; i < fNParticles; i++){

        fParticleGun->SetParticlePolarization(G4ThreeVector(2*G4UniformRand()-1,2*G4UniformRand()-1,2*G4UniformRand()-1 ) );

        //what is the particle
        fParticleGun->SetParticleDefinition(G4OpticalPhoton::OpticalPhotonDefinition());
        //determine particle momentum direction
        DirectionDecider();

        //determine particle position
        PositionDecider(currentVoxel.xPos,currentVoxel.yPos,currentVoxel.zPos,fBinWidth);
        if(fCurrentPosition == G4ThreeVector(1000000,1000000,1000000))return;// break;

        //particle direction, position, and energy sent to ParticleGun
        fParticleGun->SetParticlePosition(fCurrentPosition);
        fParticleGun->SetParticleMomentumDirection(fDirection);
        fParticleGun->SetParticleEnergy(fCurrentEnergy);
        fParticleGun->SetNumberOfParticles(1);

        //vertex generated by ParticleGun
        fParticleGun->GeneratePrimaryVertex(event);

	/*
      }
    }
  }*/
}
