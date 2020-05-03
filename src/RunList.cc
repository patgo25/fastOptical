
#include "RunList.hh"
#include "G4RunManager.hh"
#include "g4root.hh"

RunList::RunList(L200ParticleGenerator* generator, MapRunAction* mra)
	: generator(generator), mra(mra), filename("test.root")
{
 	analysis = G4Root::G4AnalysisManager::Instance();
    //openFile();

	writeDir = new G4UIdirectory("/write/");
  writeDir->SetGuidance("Output file properties");

  writeFilename= new G4UIcmdWithAString("/write/filename",this);
  writeFilename->SetGuidance("Set filename for root output file");
}


RunList::~RunList(){
	if(analysis != NULL){
		if(analysis->IsOpenFile()){
		    analysis->Write();
		    analysis->CloseFile();
		}
    	delete analysis;
	}

	//TODO: delete commands here (if any)
}


void RunList::SetNewValue(G4UIcommand *cmd, G4String newValue){
	if(cmd == writeFilename){
		filename = newValue;
	}
}

void print(MapRunAction* mra){
	for(int i = 0; i < mra->getVolumeNr(); i++){
		std::cout << "   (1) Volume "<<i<<": counts: "<<mra->getCount(i+1)<<std::endl;
	}
}

void RunList::startRuns(){
	G4RunManager* rm = G4RunManager::GetRunManager();
	openFile();
	while(true){
		G4int nrPrimaries = generator->nextVoxel();
		if(nrPrimaries == 0) break;
		rm->BeamOn(nrPrimaries);
		std::cout << " (0) run in voxel "<<generator->getCurrentVoxel()<<" ended: "<<std::endl;
		print(mra);

		writeRun(nrPrimaries);
	}
	std::cout << "Runs done "<<std::endl;
}


//only to be called ONCE
void RunList::openFile(){
    //ntuples are declared BEFORE the actual file is opened
    //see G4Simple for the reason
    //see lines 242 ff
    analysis->CreateNtuple("map","geant4 map data");
    analysis->CreateNtupleDColumn("xPos");//D for double
    analysis->CreateNtupleDColumn("yPos");
    analysis->CreateNtupleDColumn("zPos");
    analysis->CreateNtupleIColumn("counts"); //I for int
	analysis->CreateNtupleIColumn("initialNr"); //I for int
    //more if you want...

    analysis->FinishNtuple();
    analysis->SetFileName(filename);
    std::cout << "Opening file " << analysis->GetFileName() << std::endl;
    analysis->OpenFile();

    clearVars();
}

void RunList::clearVars(){
	//leave empty as we do not have vectors
}

void RunList::writeRun(G4int nrPrimaries){
	L200ParticleGenerator::Voxel voxel = generator->getCurrentVoxel();
	voxelX = voxel.xPos + 0.5*voxel.xWid;
	voxelY = voxel.yPos + 0.5*voxel.yWid;
	voxelZ = voxel.zPos + 0.5*voxel.zWid;
	count = (generator->isCurrentVoxelAborted()) ? 0 : mra->getCount(1);	//should now have only cnts in volumes with ID 1 (check macro!!!)
	initialNr = nrPrimaries;

	analysis->FillNtupleDColumn(0, voxelX);		//dont mess up ordering!
	analysis->FillNtupleDColumn(1, voxelY);
	analysis->FillNtupleDColumn(2, voxelZ);
	analysis->FillNtupleIColumn(3, count);
	analysis->FillNtupleIColumn(4, initialNr);

	analysis->AddNtupleRow();

	clearVars();
}

















