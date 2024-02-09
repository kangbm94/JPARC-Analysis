int ReadConfLine(ifstream& file, double* buf, TString& buf_l){
	if(!file.is_open()){
		cout<<"file not open"<<endl;
		return false;
	}
	TString line;
	if(file.good()&&line.ReadLine(file)){
		buf_l = line;
		if(line.IsNull()){
			return 0;
		}
		if(line[0]== '#'){
			buf[0] = -999;
			return -1;
		}
		line.ReplaceAll(",","");
		line.ReplaceAll("\"","");
		std::istringstream iss(line.Data());
		std::istream_iterator<std::string> begin(iss);
		std::istream_iterator<std::string> end;
		std::vector<TString> v(begin, end);
		int nl = v.size();
		for(int i=0;i<nl;++i){
			buf[i]=v[i].Atof();
		}
		return 1;
	}
	else{
		return 0;
	}
}
TString dir = "param/HDPRM/CH2runs/";
TString dir_org = "param/HDPRM/runs/";
TF1* fGaus = new TF1("fGaus","gaus",-1,1);
double FitHist(TH1D* hist){
	int mb = hist->GetMaximumBin();
	double peak = hist->GetXaxis()->GetBinCenter(mb);
	double w = hist->GetXaxis()->GetBinWidth(mb);
	fGaus->SetRange(peak-3*w,peak+3*w);	
	hist->Fit("fGaus","R0");
	//	cout<<"Fitted"<<endl;
	return fGaus->GetParameter(1);
}

double dat[20];
TString l;

ifstream f;
//vector<int>runlist ={5641,5642,5643,5644,5645,5646,5647,5649,5650,5652,5653,5655,5656,5657,5658,5659,5660,5661,5662,5663,5664,5665,5666}; 
vector<int>runlist ={5641,5646,5652,5659,5666}; 
void MakeBH2CableOffset(int runnum);
void CompareBH2CableOffset(int runnum);
void HodoParam(){
//	for(int run:runlist)MakeBH2CableOffset(run);
}
void MakeBH2CableOffset(int runnum){
	cout<<Form("Processing run %d",runnum)<<endl;
	TString param = Form("HodoParam_0%d",runnum);
	param = dir_org + param;
	f.open(param);
	fstream f2;
	f2.open(dir+Form("HodoParam_0%d",runnum),fstream::out);
	int flag = ReadConfLine(f,dat,l);
	vector<double> t_off;
	TFile* file = new TFile(Form("rootfiles/run0%d_Hodoscope.root",runnum));
	TCanvas* c1 = new TCanvas("c1","c1",1200,800);
	c1->Divide(4,2);
	TFile* file2 = new TFile(Form("run0%d_BH2CableHists",runnum));
	TH1D* hists[8];
	for(int i=0;i<8;++i){
		hists[i] = (TH1D*)file->Get(Form("h20%d26",i+1));
	}
	while(flag){
		if(flag<0)
			f2<<l<<endl;
		if(flag>0){
			if(dat[0]==2 and dat[1]== 0 and dat[3]==1 and dat[4]==2){//BH2ID = 2,Plane 0, seg(dat[2]), IsTDC==1, UorD == 2
					int dati[7];
					for(int i=0;i<7;++i){
						dati[i]=dat[i];
					}
					int seg = dati[2];
					cout<<seg<<endl;
					c1->cd(seg+1);
					double t_offset = FitHist(hists[seg]);
					t_offset+=dat[5];
					TString param = Form("%d %d %d %d %d %f %d",dati[0],dati[1],dati[2],dati[3],dati[4],t_offset,dati[6]);
					f2<<param.Data()<<endl;
			}
			else 
				f2<<l<<endl;
		}
		flag = ReadConfLine	(f,dat,l);
	}
	f2.close();
	f.close();
	delete file;
	delete c1;
}
void CompareBH2CableOffset(){
	TString rdir = dir;
	TFile* file = new TFile("BH2Offset.root","recreate");
	TTree* tree = new TTree("tree","tree");
	int rn = 0;
	double t_off[8];
	tree->Branch("runnum",&rn);
	tree->Branch("BH2Offset",t_off,"BH2Offset[8]/D");
	for(auto run:runlist){
		TString param = Form("HodoParam_0%d",run);
		param = rdir + param;
		ifstream pfile;
		pfile.open(param);
		int flag = ReadConfLine(pfile,dat,l);
		rn = run;
		while(flag){
			if(flag>0){
				if(dat[0]==2 and dat[1]== 0 and dat[3]==1 and dat[4]==2){//BH2ID = 2,Plane 0, seg(dat[2]), IsTDC==1, UorD == 2
					int dati[7];
					for(int i=0;i<7;++i){
						dati[i]=dat[i];
					}
					int seg = dati[2];
					cout<<seg<<endl;
					double t_offset=dat[5];
					t_off[seg]=t_offset;
				}
			}
		flag = ReadConfLine(pfile,dat,l);
		}
		tree->Fill();
	}
	file->Write();
}
