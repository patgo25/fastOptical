#ifndef L200ParticleGeneratorMessenger_hh
#define L200ParticleGeneratorMessenger_hh

//---------------------------------------------------------------------------//

#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"

//---------------------------------------------------------------------------//

class G4UIdirectory;
class G4UIcmdWithAString;
class l200ParticleGenerator;

class L200ParticleGeneratorMessenger : public G4UImessenger
{
public:
	L200ParticleGeneratorMessenger(L200ParticleGenerator *generator);
	~L200ParticleGeneratorMessenger();

	void SetNewValue(G4UIcommand *cmd, G4String newValue);

protected:

private:
  L200ParticleGenerator* fLiquidArgonGenerator;
  G4UIdirectory* fLiquidArgonDirectory;
  G4UIcmdWithADoubleAndUnit* fLiquidArgonSetRadius;
  G4UIcmdWithADoubleAndUnit* fLiquidArgonSetRadiusMin;
  G4UIcmdWithADoubleAndUnit* fLiquidArgonSetHeight;
  G4UIcmdWith3VectorAndUnit* fLiquidArgonSetCenterVector;
  G4UIcmdWithADoubleAndUnit* fLiquidArgonSetBinWidth;
  G4UIcmdWithADouble* fLiquidArgonSetNParticles;
  G4UIcmdWithABool* fLiquidArgonSet1D;

};
#endif
