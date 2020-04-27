/*
* Root script, that takes raw counts-in-voxels data from L200 optical simulation and 
* creates fancy 2D histos and stuff.
*
*
*
*/

TH2D* getSegment(){
	TTree* map = (TTree*) gFile->Get("map");
	//map->Show(0);

	double xPos, yPos, zPos;
	int counts, initialNr;

	map->SetBranchAddress("xPos", &xPos);
	map->SetBranchAddress("yPos", &yPos);
	map->SetBranchAddress("zPos", &zPos);
	map->SetBranchAddress("counts", &counts);
	map->SetBranchAddress("initialNr", &initialNr);
	
	double xMin = 1e50, xMax = -1e50, yMin = 1e50, yMax = -1e50;
	double xPitch = 1e50, yPitch = 1e50;

	int entries = map->GetEntries();
	for(int i = 0; i < entries; i++){
		map->GetEntry(i);

		if(xPos > xMax) xMax = xPos;
		if(xPos < xMin) xMin = xPos;
		if(yPos > yMax) yMax = yPos;
		if(yPos < yMin) yMin = yPos;
	}	

	//search for points with minimum (but still non-zero) distance to leftmost point
	// --> gives pitch (const distance btw 2 points)
	for(int i = 0; i < entries; i++){
		map->GetEntry(i);

		double xDiff = xPos - xMin;	//>= 0 by construction
		double yDiff = yPos - yMin; //>= 0 by construction

		if(xDiff > 1e-6 && xDiff < xPitch) xPitch = xDiff;
		if(yDiff > 1e-6 && yDiff < yPitch) yPitch = yDiff;
	}	

	std::cout << xMin<<", "<<yMin<<", "<<xPitch<<", "<<yPitch<<std::endl;

	int xBins = (int)((xMax-xMin)/xPitch+0.5);
	int yBins = (int)((yMax-yMin)/yPitch+0.5);

	TH2D* histo = new TH2D("segment", "2D map segment", xBins, xMin, xMax, yBins, yMin, yMax);

	for(int i = 0; i < entries; i++){
		map->GetEntry(i);

		histo->SetBinContent(histo->FindBin(xPos,yPos),((double)counts)/initialNr);
	}	
	

	return histo;
}

























