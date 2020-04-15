#indef SENSITIVE_DETECTOR_H
#define SENSITIVE_DETECTOR_H

#include "G4VSensitiveDetector.hh"  
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4UnitsTable.hh"

#include "Hit.hh"
#include <time.h>

#include <fstream>

class MySensitiveDetector : public G4VSensitiveDetector
{
  public:
      MySensitiveDetector(G4String);
     ~MySensitiveDetector();

      //! create an instance of ExN06LDHitsCollection and add a new hits collection
      //! to the G4HCofThisEvent instance
      void Initialize(G4HCofThisEvent* hitCollection);

      G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);


	//! write the LD Data to the root file if necessary
	void EndOfEvent(G4HCofThisEvent* hitCollection);

      //! Return the number of Hits
      G4int getEntries();
	G4int numbersOfPhotons;
	G4double PE;
		/*

  private:
	MyHitsCollection* HitsCollection;
	std::ofstream (*ResultFile);
	std::ofstream (*ResultFile2);
 
 protected:
	G4CollectionNameVector collectionName;*/


};

#endif
