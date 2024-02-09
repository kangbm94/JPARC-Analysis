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
TF1* fGaus = new TF1("fGaus","gaus",-1,1);
double FitHist(TH1D* hist){
	int mb = hist->GetMaximumBin();
	double peak = hist->GetXaxis()->GetBinCenter(mb);
	double w = hist->GetXaxis()->GetBinWidth(mb);
	fGaus->SetRange(peak-3*w,peak+3*w);	
	hist->Fit("fGaus","R");
	//	cout<<"Fitted"<<endl;
	return fGaus->GetParameter(1);
}

double dat[20];
TString l;

string param = "HodoParam_KBM";
ifstream f;
void HodoParam(){
	f.open(param);
}
void MakeBH2CableOffset(int runnum){
	fstream f2;
	f2.open(Form("HodoParam_0%d",runnum),fstream::out);
	int flag = ReadConfLine(f,dat,l);
	vector<double> t_off;
	TH1D* hists[8];
	TFile* file = new TFile(Form("run0%d_Hodoscope.root",runnum));
	TCanvas* c1 = new TCanvas("c1","c1",1200,800);
	c1->Divide(4,2);
	for(int i=0;i<8;++i){
		hists[i] = (TH1D*)file->Get(Form("h20%d19",i+1));
	}
	while(flag){
		if(flag<0)
			f2<<l<<endl;
		if(flag>0){
			if(dat[0]==2 and dat[1]== 0 and dat[3]==1 and dat[4]==2){//BH2ID = 2,Plane 0, seg(dat[2]), IsTDC==1, UorD == 2
					t_off.push_back(dat[5]);
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
	for(auto b : t_off) cout<<b<<" ns"<<endl;
}

