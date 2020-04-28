
#ifndef L200ParticleGenerator_h
#define L200ParticleGenerator_h

//---------------------------------------------------------------------------//

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"


//---------------------------------------------------------------------------//

class G4Event;
class L200ParticleGeneratorMessenger;
class G4ParticleGun;
class G4Run;

class L200ParticleGenerator
{

  public:

	struct Voxel{
		G4double xPos, yPos, zPos;	//all lower points; i.e. x from xPos to xWid, ...
		G4double xWid, yWid, zWid;		//TODO currently unused; need to be implemented for oblong vxels

		friend std::ostream& operator<<(std::ostream& os, const Voxel& vx){
			os << "("<<vx.xPos<<","<<vx.yPos<<","<<vx.zPos<<")";
			return os;
		};
	};

    L200ParticleGenerator();
    ~L200ParticleGenerator();

    void GeneratePrimaryVertex(G4Event *event);
    void SetParticlePosition(G4ThreeVector pos) { fCurrentPosition = pos;}
    void DirectionDecider();
    void PositionDecider(G4double x,G4double y,G4double z,G4double binWidth);
    G4bool IsInArgon(G4ThreeVector rp);

    //Messenger Commands
    void SetRadius(G4double r){fRadiusMax = r;}
    void SetRadiusMin(G4double r){fRadiusMin = r;}
    void SetHeight(G4double h){fZ = h;}
    void SetCenterVector(G4ThreeVector vec){fCenterVector = vec;}
    void SetBinWidth(G4double width) {fBinWidth = width;}
    void SetNParticles(G4double N) {fNParticles = N;}
    void is1DScan(G4bool b){fis1D =b;}
	void setVerbosity(G4int verbose){verbosity = verbose;};
	void setAbortOnNonlar(G4bool flag){abortOnNonlar = flag;};

	//methods called from above (RunList)
	//to be called BEFORE STARTING ANY RUN!!!
	int nextVoxel();		//commands the generator to switch to the next voxel (set internals accordingly)
							//return is nr of primaries to be shot in there (beamOn called outside).
							// return 0 --> end runs

	Voxel getCurrentVoxel() {return currentVoxel;};
	G4bool isCurrentVoxelAborted(){return abortVoxel;};

  private:
    static const G4double LambdaE;
    G4ParticleGun*	fParticleGun;
	G4int verbosity;	
	
    //particle properties
    G4double  fCurrentEnergy; // energy of current particle
    G4ThreeVector fCurrentPosition; // current position of particle
    G4ThreeVector fDirection; // direction of momentum
    G4double fRadiusMax = 0;		//radius directly affects x from/to; y borders 
								//dynamically calculated from xPos & scanAngle
    G4double fRadiusMin = 0;
    G4double fZ = 0;
    G4double fBinWidth = 0;
	G4double scanAngle;		//defines y borders as from 0 to at least x*tan(scanAngle)
		
    G4double fNParticles = 1;
    G4ThreeVector fCenterVector;
    G4String fParticleType = "opticalphoton";
    G4double fEnergy = 0;
    G4bool fis1D = true;
    G4double fLArWL;
	G4bool abortOnNonlar;		//if a single non-LAr hit should abort the whole voxel

    L200ParticleGeneratorMessenger* fMessenger;

	Voxel currentVoxel;
	G4bool larFailed;		//set in PositionDecider when no proper LAr position found
	G4bool abortVoxel;		//set in Generate... @ hit of nonLar, read in Generate... & from RunList,
							//reset in nextVoxel

	uint32_t flatVoxelIndex;	//flat index for current voxel (defines both x and y index if needed)

};
#endif
