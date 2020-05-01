/*
* Root script, that takes raw counts-in-voxels data from L200 optical simulation and
* creates fancy 2D histos and stuff.
*
*
*
*/
void pad2png(TH2D* h,bool is1D=false)
{
	TCanvas* c = new TCanvas();
	c->SetWindowSize(600*2,150*2);
	if(is1D){

	}
   	h->Draw("colz");
	h->SetMaximum(0.04);
	gStyle->SetOptStat(0);
	string title = h->GetTitle();
	title.append(".png");
	gSystem->ProcessEvents();

 	TImage *img = TImage::Create();
	//img->FromPad(c, 10, 10, 300, 200);
 	img->FromPad(c);
	img->WriteImage(title.c_str());

	delete c;
	delete img;
}


TH2D* getSegment(const char* file){
	TFile* theFile = new TFile(file);
	TTree* map = (TTree*) theFile->Get("map");
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
	cout << xBins << endl;
	TH2D* histo = new TH2D("segment", "2D map segment", xBins, xMin, xMax, yBins, yMin, yMax);

	for(int i = 0; i < entries; i++){
		map->GetEntry(i);

		histo->SetBinContent(histo->FindBin(xPos,yPos),((double)counts)/initialNr);
	}


	return histo;
}

TH2D* getRotSymPlot(const char* file, double cuts=14 ){

	TH2D* cakepiece = NULL;
	cakepiece = getSegment(file);

	double binwidth = cakepiece->GetXaxis()->GetBinWidth(1);
	double xCakepieceMax = cakepiece->GetXaxis()->GetXmax();
	int xBins = cakepiece->GetNbinsX();
	int yBins = cakepiece->GetNbinsY();
	cout << "num of bins: x= " << xBins <<"XMax: "<< xCakepieceMax << " y= " << yBins <<" width= " << binwidth << " angle cuts = " << 360./cuts << endl;
	double cakeMin = -1*xCakepieceMax;
	int cakeBins = (int)(xCakepieceMax-cakeMin)/binwidth;

	TH2D* cake = new TH2D("cake","cake",cakeBins, cakeMin,xCakepieceMax,cakeBins,cakeMin,xCakepieceMax);
	for(int i=0; i<=cakeBins;i++){
		double x = cake->GetXaxis()->GetBinCenter(i);
		for(int j=0;j<=cakeBins;j++){
			double y = cake->GetYaxis()->GetBinLowEdge(j);
			if(x*x +y*y > xCakepieceMax*xCakepieceMax) continue;
			double counts = 0;
			for(int angi=0; angi<cuts;angi++){
				double ang =angi*(2.*TMath::Pi()/cuts);
				double rotx = x*cos(ang)-y*sin(ang);
				double roty = x*sin(ang) + y*cos(ang);
				counts = cakepiece->GetBinContent(cakepiece->FindBin(rotx,roty));
				if(counts>0) break;
			}
			if(counts == 0){
				for(int angi=0; angi<cuts;angi++){
					double ang =angi*(2.*TMath::Pi()/cuts);
					double rotx = x*cos(ang)-y*sin(ang);
					double roty = x*sin(ang) + y*cos(ang);
					counts = cakepiece->GetBinContent(cakepiece->FindBin(rotx,-roty));
					if(counts>0) break;
				}

			}
			cake->SetBinContent(i,j,counts);
		}
	}

	return cake;
}

void makeGifPlot(const char* dirname, bool is1D=false,bool cake=false){

	TSystemDirectory dir(dirname, dirname);
        TList *files = dir.GetListOfFiles();
        if(files){
                TSystemFile *file;
                TString fname;
                TIter next(files);
		TH2D* plot = NULL;
                while((file=(TSystemFile*)next())){
                        fname = file->GetName();
                        if(!file->IsDirectory() && fname.EndsWith(".root")){
				string path = dirname;
				path.append("/");
				path.append(fname);
                                cout << path<< endl;
                                cout <<"Processing Absorption length of ";
                                TString lengths = fname(fname.Last('-')+1,fname.First('.')-fname.Last('-')-1);
                                int length = std::stoi(lengths.Data());
				cout << length << endl;
				if(cake)
					plot=getRotSymPlot(path.c_str());
				else
					plot = getSegment(path.c_str());
				if(length < 1000)
					{string title = "0"; title.append(lengths); lengths=title;}
				plot->SetTitle(lengths);
				pad2png(plot);
                        }
                }
        }
}

TGraph* get1DPlot(const char* file, const char* title){

	TFile fileIn(file);
  	TTree* map = nullptr;
   	fileIn.GetObject("map",map);
	int n = map->Draw("xPos:counts/initialNr","","goff");
	TGraph *g = new TGraph(n,map->GetV1(),map->GetV2());
	g->GetXaxis()->SetTitle("Radius [mm]");
	g->GetYaxis()->SetTitle("Detection probability");
	g->SetTitle(title);
	return g;
}


























