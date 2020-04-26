#ifndef L200BoundaryProcess_h
#define L200BoundaryProcess_h


// Includes

#include "globals.hh"
#include "templates.hh"
#include "geomdefs.hh"
#include "Randomize.hh"

#include "G4RandomTools.hh"
#include "G4RandomDirection.hh"

#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4OpticalPhoton.hh"
#include "G4TransportationManager.hh"

// Class Description:
// Discrete Process -- reflection/refraction at optical interfaces.
// Class inherits publicly from G4VDiscreteProcess.
// Class Description - End:

// Class Definition
enum L200OpBoundaryProcessStatus {  Undefined,
                                  FresnelRefraction, FresnelReflection,
                                  TotalInternalReflection,
                                  LambertianReflection, LobeReflection,
                                  SpikeReflection, BackScattering,
                                  Absorption, Detection, NotAtBoundary,
                                  SameMaterial, StepTooSmall, NoRINDEX,
                                  PolishedLumirrorAirReflection,
                                  PolishedLumirrorGlueReflection,
                                  PolishedAirReflection,
                                  PolishedTeflonAirReflection,
                                  PolishedTiOAirReflection,
                                  PolishedTyvekAirReflection,
                                  PolishedVM2000AirReflection,
                                  PolishedVM2000GlueReflection,
                                  EtchedLumirrorAirReflection,
                                  EtchedLumirrorGlueReflection,
                                  EtchedAirReflection,
                                  EtchedTeflonAirReflection,
                                  EtchedTiOAirReflection,
                                  EtchedTyvekAirReflection,
                                  EtchedVM2000AirReflection,
                                  EtchedVM2000GlueReflection,
                                  GroundLumirrorAirReflection,
                                  GroundLumirrorGlueReflection,
                                  GroundAirReflection,
                                  GroundTeflonAirReflection,
                                  GroundTiOAirReflection,
                                  GroundTyvekAirReflection,
                                  GroundVM2000AirReflection,
                                  GroundVM2000GlueReflection,
				  TPBMagic };

class L200OpBoundaryProcess : public G4VDiscreteProcess
{

public:

        // Constructors and Destructor

        L200OpBoundaryProcess(const G4String& processName = "OpBoundary",
                                     G4ProcessType type = fOptical);
        ~L200OpBoundaryProcess();

private:

        L200OpBoundaryProcess(const L200OpBoundaryProcess &right);

        // Operators

        L200OpBoundaryProcess& operator=(const L200OpBoundaryProcess &right);

public:

        // Methods

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable' only for an optical photon.

        G4double GetMeanFreePath(const G4Track& ,
                                 G4double ,
                                 G4ForceCondition* condition);
        // Returns infinity; i. e. the process does not limit the step,
        // but sets the 'Forced' condition for the DoIt to be invoked at
        // every step. However, only at a boundary will any action be
        // taken.

        G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                       const G4Step&  aStep);
        // This is the method implementing boundary processes.

        L200OpBoundaryProcessStatus GetStatus() const;
        // Returns the current status.

	//Magic setter
	void setFiberHitProb(G4double value){theProb = value;}
	void setMagicMaterialName(G4String value){theTPBMagicMaterialName = value;}
	void setLArWL(G4double value){theLArWL = value;}

private:

        G4bool G4BooleanRand(const G4double prob) const;

        G4ThreeVector GetFacetNormal(const G4ThreeVector& Momentum,
                                     const G4ThreeVector&  Normal) const;

        void DielectricMetal();
        void DielectricDielectric();
        void DielectricLUT();

        void ChooseReflection();
        void DoAbsorption();
        void DoReflection();

        G4double GetIncidentAngle();
        // Returns the incident angle of optical photon

        G4double GetReflectivity(G4double E1_perp,
                                 G4double E1_parl,
                                 G4double incidentangle,
                                 G4double RealRindex,
                                 G4double ImaginaryRindex);
        // Returns the Reflectivity on a metalic surface

        void CalculateReflectivity(void);

        void BoundaryProcessVerbose(void) const;

private:

        G4double thePhotonMomentum;

        G4ThreeVector OldMomentum;
        G4ThreeVector OldPolarization;

        G4ThreeVector NewMomentum;
        G4ThreeVector NewPolarization;

        G4ThreeVector theGlobalNormal;
        G4ThreeVector theFacetNormal;

        G4Material* Material1;
        G4Material* Material2;

        G4OpticalSurface* OpticalSurface;

        G4MaterialPropertyVector* PropertyPointer;
        G4MaterialPropertyVector* PropertyPointer1;
        G4MaterialPropertyVector* PropertyPointer2;

        G4double Rindex1;
        G4double Rindex2;

        G4double cost1, cost2, sint1, sint2;

        L200OpBoundaryProcessStatus theStatus;

        G4OpticalSurfaceModel theModel;

        G4OpticalSurfaceFinish theFinish;

        G4double theReflectivity;
        G4double theEfficiency;
        G4double theTransmittance;
        G4double prob_sl, prob_ss, prob_bs;

        G4int iTE, iTM;

        //Magic TPB variables
	G4double kCarTolerance;
	G4double theProb;
	G4String theTPBMagicMaterialName;
	G4double theLArWL;

};

// Inline methods

inline
G4bool L200OpBoundaryProcess::G4BooleanRand(const G4double prob) const
{
  /* Returns a random boolean variable with the specified probability */

  return (G4UniformRand() < prob);
}

inline
G4bool L200OpBoundaryProcess::IsApplicable(const G4ParticleDefinition&
                                                       aParticleType)
{
   return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
}

inline
L200OpBoundaryProcessStatus L200OpBoundaryProcess::GetStatus() const
{
   return theStatus;
}

inline
void L200OpBoundaryProcess::ChooseReflection()
{
                 G4double rand = G4UniformRand();
                 if ( rand >= 0.0 && rand < prob_ss ) {
                    theStatus = SpikeReflection;
                    theFacetNormal = theGlobalNormal;
                 }
                 else if ( rand >= prob_ss &&
                           rand <= prob_ss+prob_sl) {
                    theStatus = LobeReflection;
                 }
                 else if ( rand > prob_ss+prob_sl &&
                           rand < prob_ss+prob_sl+prob_bs ) {
                    theStatus = BackScattering;
                 }
                 else {
                    theStatus = LambertianReflection;
                 }
}

inline
void L200OpBoundaryProcess::DoAbsorption()
{
              theStatus = Absorption;

              if ( G4BooleanRand(theEfficiency) ) {

                 // EnergyDeposited =/= 0 means: photon has been detected
                 theStatus = Detection;
                 aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
              }
              else {
                 aParticleChange.ProposeLocalEnergyDeposit(0.0);
              }

              NewMomentum = OldMomentum;
              NewPolarization = OldPolarization;

//              aParticleChange.ProposeEnergy(0.0);
              aParticleChange.ProposeTrackStatus(fStopAndKill);
}

inline
void L200OpBoundaryProcess::DoReflection()
{
        if ( theStatus == LambertianReflection ) {

          NewMomentum = G4LambertianRand(theGlobalNormal);
          theFacetNormal = (NewMomentum - OldMomentum).unit();

        }
        else if ( theFinish == ground ) {

          theStatus = LobeReflection;
          if ( PropertyPointer1 && PropertyPointer2 ){
          } else {
             theFacetNormal =
                 GetFacetNormal(OldMomentum,theGlobalNormal);
          }
          G4double PdotN = OldMomentum * theFacetNormal;
          NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;

        }
        else {

          theStatus = SpikeReflection;
          theFacetNormal = theGlobalNormal;
          G4double PdotN = OldMomentum * theFacetNormal;
          NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;

        }
        G4double EdotN = OldPolarization * theFacetNormal;
        NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;
}

#endif
