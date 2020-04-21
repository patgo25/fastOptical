
#include "RunList.hh"
#include "G4RunManager.hh"

RunList::RunList(L200ParticleGenerator* generator, MapRunAction* mra)
	: generator(generator), mra(mra)
{

}


RunList::~RunList(){
	//TODO: delete commands here (if any)
}


void RunList::SetNewValue(G4UIcommand *cmd, G4String newValue){}

void print(MapRunAction* mra){
	for(int i = 0; i < mra->getVolumeNr(); i++){
		std::cout << "   (1) Volume "<<i<<": counts: "<<mra->getCount(i+1)<<std::endl;
	}
}

void RunList::startRuns(){
	G4RunManager* rm = G4RunManager::GetRunManager();
	while(true){
		int nrPrimaries = generator->nextVoxel();
		if(nrPrimaries == 0) break;
		rm->BeamOn(nrPrimaries);
		std::cout << " (0) run in voxel "<<generator->getCurrentVoxel()<<" ended: "<<std::endl;
		print(mra);
	}
	std::cout << "Runs done "<<std::endl;
}


