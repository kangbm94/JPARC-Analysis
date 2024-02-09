#include "TPCPadHelper.hh"
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
TF1* fGaus = new TF1("fGaus","gaus",-15,15);
double FitHist(TH1D* hist,double sigw = 3){
	int mb = hist->GetMaximumBin();
	double peak = hist->GetXaxis()->GetBinCenter(mb);
	double w = hist->GetXaxis()->GetBinWidth(mb);
	fGaus->SetRange(peak-sigw*w,peak+sigw*w);	
	hist->Fit("fGaus","R");
	//	cout<<"Fitted"<<endl;
	return fGaus->GetParameter(1);
}

double dat[20];
TString line;

string param = "TPCPositionCorrectionMap_230310";
ifstream f;
int runnum = 5856;
void TPCPos(){
	f.open(param);
	cout<<"CorrectTPCPos()"<<endl;
}
void CorrectTPCPos(){
	gStyle->SetOptFit(10);
	fstream f2;
	f2.open(Form("TPCPositionCorrectionMap_KBM"),fstream::out);
	int flag = ReadConfLine(f,dat,line);
	TFile* file = new TFile(Form("NoCorL29/run0%d_Beamthrough.root",runnum));
//	TFile* file = new TFile(Form("AllCor/run0%d_Beamthrough.root",runnum));
//	TFile* file = new TFile(Form("NoCorL29/Beamthrough.root"));
	TTree* tree = (TTree*)file->Get("tpc");
	double mom0;
	int Ylist[5] = {-250,-150,-50,50,150};
	vector<double>*layer = new vector<double>;
	vector<double>*row = new vector<double>;
	vector<double>*hitpos_x = new vector<double>;
	vector<double>*hitpos_y = new vector<double>;
	vector<double>*hitpos_z = new vector<double>;
	vector<double>*residual_x = new vector<double>;
	vector<double>*residual_y = new vector<double>;
	vector<double>*residual_z = new vector<double>;
	row = new vector<double>;
	tree->SetBranchAddress("mom0",&mom0);
	tree->SetBranchAddress("layer",&layer);
	tree->SetBranchAddress("row",&row);
	tree->SetBranchAddress("hitpos_x",&hitpos_x);
	tree->SetBranchAddress("hitpos_y",&hitpos_y);
	tree->SetBranchAddress("hitpos_z",&hitpos_z);
	tree->SetBranchAddress("residual_x",&residual_x);
	tree->SetBranchAddress("residual_y",&residual_y);
	tree->SetBranchAddress("residual_z",&residual_z);
//	TFile* file2 = new TFile(Form("run0%d_BeamthroughHists.root",runnum),"recreate");
	TFile* file2 = new TFile(Form("BeamthroughHists.root"),"recreate");
	TTree* tree2 = new TTree("tree","tree");
	TH1D* ResHistX[32][244];
	TH1D* ResHistY[32][244];
	TH1D* ResHistZ[32][244];
	int corL,corR;
	double dX,dZ,rX,rZ;
	tree2->Branch("layer",&corL);
	tree2->Branch("row",&corR);
	tree2->Branch("dX",&dX);
	tree2->Branch("dZ",&dZ);
	tree2->Branch("rX",&rX);
	tree2->Branch("rZ",&rZ);
	for(int layer =0;layer<32;++layer){
		int max_row = tpc::padParameter[layer][1];
		for(int row =0;row<max_row;++row){
			TString titlex = Form("ResXL%d_R%d",layer,row);
			TString titley = Form("ResYL%d_R%d",layer,row);
			TString titlez = Form("ResZL%d_R%d",layer,row);
			ResHistX[layer][row] = new TH1D(titlex,titlex,150,-15,15);
			ResHistY[layer][row] = new TH1D(titley,titley,150,-15,15);
			ResHistZ[layer][row] = new TH1D(titlez,titlez,150,-15,15);
//			cout<<titlex<<endl;
//			cout<<ResHistX[layer][row]->GetEffectiveEntries()<<endl;
		}
	}
	int ent = tree->GetEntries();
	for(int iev=0;iev<ent;++iev){
		tree->GetEntry(iev);
		int nh = layer->size();
		for(int ih=0;ih<nh;++ih){
			int l =(int) layer->at(ih);
			int r =(int) row->at(ih);
			double rx = residual_x->at(ih);
			double ry = residual_y->at(ih);
			double rz = residual_z->at(ih);
//			cout<<l<<" , "<<r<<endl;
//			cout<<rx<<" , "<<ry<<" , "<<rz<<endl;
//			TString titlex = Form("ResXL%d_R%d",l,r);
//			cout<<titlex<<endl;
			ResHistX[l][r]->Fill(rx);
			ResHistY[l][r]->Fill(ry);
			ResHistZ[l][r]->Fill(rz);
		}
	}
	int nl = 0;
	while(flag){
		if(nl%1000==0) cout<<"Line "<<nl<<endl;
		nl++;
		if(flag<0) f2<<line<<endl;
		if(flag>0){
			if(dat[0] < 29){
				f2<<line<<endl;
			}
			else{
				int lay = dat[0];
//				int max_row = tpc::padParameter[lay][1];
				int rw = dat[1];
				if(ResHistX[lay][rw]->GetEffectiveEntries()<20){
					f2<<line<<endl;
					flag = ReadConfLine	(f,dat,line);
					continue;
				}
				double Mean = ResHistX[lay][rw]->GetMean();
				double pb = ResHistX[lay][rw]->GetMaximumBin();
				if(ResHistX[lay][rw]->GetMaximum()<8){
					f2<<line<<endl;
					flag = ReadConfLine	(f,dat,line);
					continue;
				}

				double Peak=ResHistX[lay][rw]->GetBinCenter(pb);
				double W = ResHistX[lay][rw]->GetBinWidth(50);
				fGaus->SetRange(10,10);
//				fGaus->SetParLimits(1,Peak-3*W,Peak+3*W);
				fGaus->SetParLimits(1,Peak - 2*W,Peak + 2*W);
				fGaus->SetParLimits(2,0.1,0.6);
//				ResHistX[lay][rw]->Fit("fGaus","R");
//				Peak = fGaus->GetParameter(1);
				fGaus->SetRange(Peak-3*W,Peak+3*W);
				ResHistX[lay][rw]->Fit("fGaus","R");
				rX = fGaus->GetParameter(1);
				Mean = ResHistZ[lay][rw]->GetMean();
				W = ResHistZ[lay][rw]->GetBinWidth(50);
				fGaus->SetRange(10,10);
//				fGaus->SetRange(Mean-3*W,Mean+3*W);
				fGaus->SetParLimits(1,-2,2);
				fGaus->SetParLimits(2,0.01,1);
				ResHistZ[lay][rw]->Fit("fGaus","R");
				Peak = fGaus->GetParameter(1);
				fGaus->SetRange(Peak-20*W,Peak+20*W);
				ResHistZ[lay][rw]->Fit("fGaus","R");
				rZ = fGaus->GetParameter(1);
				int Y = dat[2];
				dX = dat[3]-rX; 
				dZ = dat[4]-rZ;
				//						if(abs(rX - dX) > 5 or abs(rZ - dZ)>1) continue;
				corL = lay;
				corR = rw;
				tree2->Fill();
				TString par = Form("%d %d %d %f %f",lay,rw,Y,rX,rZ);
				f2<<par<<endl;
			}
		}
		flag = ReadConfLine	(f,dat,line);
	}
	file2->Write();
}

