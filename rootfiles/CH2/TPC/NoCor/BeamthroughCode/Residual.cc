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


TF1* ResDe = new TF1("ResDe","[0]+[1]/x",90,3000);
TF1* fGaus = new TF1("fGaus","gaus",-10,10);
TF1* ResDl = new TF1("ResDl","[0]+[1]*abs(x)",-300,300);
TF1* ResDl2 = new TF1("ResDl2","[0]+[1]*x*x",-300,300);
//ToDoList : CheckLinearTracking in 3D display
//int runnum = 5721;
//int runnum = 5764;
//int runnum = 5641;
//int runnum = 5844;//5847,5855,5856,5858,5860,5864,5866
int runnum = 5721;
//int runnum = 5855;
vector<double> SmoothOut(vector<double>resi,double th){
	vector<double>ret;
	ret.push_back(resi.at(0));
	int nv = resi.size();
	for(int i = 1; i< nv-1; ++i){
		double prev = resi.at(i-1);
		double now = resi.at(i);
		double after = resi.at(i+1);
		if(prev*now*after==0){
			if(prev*after!=0 and abs(prev - after)< 1.5*th){
				ret.push_back((prev+after)/2);
			}
			else{
				cout<<Form("No Continuation... %g,%g,%g",prev,now,after)<<endl;
				ret.push_back(now);
			}
			continue;
		}
		if(abs(prev-now)<th or abs(after-now)<th){
			ret.push_back(now);
			continue;
		}
		ret.push_back((prev+after)/2);
		ret.push_back(after);
		i++;
	}
	ret.push_back(resi.at(nv-1));
	return ret;
}
vector<double> SmoothIn(vector<double>resi,double th){
	vector<double>ret;
	int nv = resi.size();
	for(int i = 0;i<nv;++i){
		double prev,after;
		if(i==0){
			prev = resi.at(nv-1);
		}
		else{
			prev = resi.at(i-1);
		}
		if(i==nv-1){
			after = resi.at(0);
		}
		else{
			after = resi.at(i+1);
		}
		double now = resi.at(i);
		if(prev*now*after==0){
			if(prev*after!=0 and abs(prev - after)< 1.5*th) ret.push_back((prev+after)/2);
			else ret.push_back(now);
			continue;
		}
		if(abs(prev-now)<th or abs(after-now)<th){
			ret.push_back(now);
			continue;
		}
		ret.push_back((prev+after)/2);
		if(i!=nv-1){
			ret.push_back(after);
			i++;
		}
	}
	return ret;
}
vector<double>Smooth(vector<double>resi,double th,bool Inner){
	if(Inner) return SmoothIn(resi,th);
	else return SmoothOut(resi,th);
}
void MakePadResidualParameter(){
//	TFile* file = new TFile(Form("run0%d_ResidualHists.root",runnum));
//	TFile* file2 = new TFile(Form("run0%d_ResidualHistsFit.root",runnum),"recreate");
	TFile* file = new TFile("BeamthroughHists.root");
	TFile* file2 = new TFile("BeamthroughHistsFit.root","recreate");
	TTree* tree = new TTree("tree","tree");
	TH1D* histX[32][244];
	TH1D* histY[32][244];
	TH1D* histZ[32][244];
	TH1D* histR[32][244];
	TH1D* histL[32][244];
	TH1D* histW[32][244];
	double peak,width;
	double rx,ry,rz,rw,rl,rxCor,rzCor,rwCor,rr,sx,sy,sz,sw,sl,sr,angle;
	double CorX,CorZ;
	int l,r,ent;
	bool isDead;
	tree->Branch("layer",&l);
	tree->Branch("row",&r);
	tree->Branch("angle",&angle);
	tree->Branch("ent",&ent);
	tree->Branch("residual_x",&rx);
	tree->Branch("residual_y",&ry);
	tree->Branch("residual_z",&rz);
	tree->Branch("residual_r",&rr);
	tree->Branch("residual_w",&rw);
	tree->Branch("residual_l",&rl);
	tree->Branch("residual_xCor",&rxCor);
	tree->Branch("residual_zCor",&rzCor);
	tree->Branch("residual_wCor",&rwCor);
	tree->Branch("resolution_x",&sx);
	tree->Branch("resolution_y",&sy);
	tree->Branch("resolution_z",&sz);
	tree->Branch("resolution_r",&sr);
	tree->Branch("resolution_w",&sw);
	tree->Branch("resolution_l",&sl);
	tree->Branch("Correction_x",&CorX);
	tree->Branch("Correction_z",&CorZ);
	tree->Branch("isDead",&isDead);

	fstream f2;
	f2.open(Form("TPCPositionCorrectionMap_KBM"),fstream::out);
	f2<<Form("%d	%d	%d",1,-250,500)<<endl;
	vector<vector<double>> AllResX;
	vector<vector<double>> AllResY;
	vector<vector<double>> AllResZ;
	vector<vector<double>> AllResR;
	vector<vector<double>> AllResW;
	vector<vector<double>> AllResL;
	vector<vector<double>> AllResXCor;
	vector<vector<double>> AllResYCor;
	vector<vector<double>> AllResZCor;
	vector<vector<double>> AllResWCor;
	vector<vector<double>> AllSigX;
	vector<vector<double>> AllSigY;
	vector<vector<double>> AllSigZ;
	vector<vector<double>> AllSigR;
	vector<vector<double>> AllSigW;
	vector<vector<double>> AllSigL;
	vector<vector<double>> ParX;
	vector<vector<double>> ParZ;
	ifstream fpar;
	fpar.open("TPCPositionCorrectionMap_KBM1st");
	TString bufl;
	double buf[10];
	int read_flag = ReadConfLine(fpar,buf,bufl);
	ParX.resize(32);
	ParZ.resize(32);
	for(int layer=0;layer<32;++layer){
		int max_row = tpc::padParameter[layer][1];
		ParX[layer].resize(max_row);
		ParZ[layer].resize(max_row);
	}

	while(read_flag){
		if(buf[1] <0){
			read_flag = ReadConfLine(fpar,buf,bufl);
			continue;
		}
		else{
			int layer = buf[0];
			int row = buf[1];
			double px  = buf[3];
			double pz  = buf[4];
			ParX[layer][row] = px;
			ParZ[layer][row] = pz;
			read_flag = ReadConfLine(fpar,buf,bufl);
		}
	}
	for(int layer =0;layer < 32;++layer){
		int max_row = tpc::padParameter[layer][1];
		vector<double> res_x;
		vector<double> res_y;
		vector<double> res_z;
		vector<double> res_r;
		vector<double> res_w;
		vector<double> res_l;
		vector<double> sig_x;
		vector<double> sig_y;
		vector<double> sig_z;
		vector<double> sig_r;
		vector<double> sig_w;
		vector<double> sig_l;
		for(int row = 0;row < max_row; ++row){
			TString titleX = Form("ResXL%d_R%d",layer,row);
			TString titleY = Form("ResYL%d_R%d",layer,row);
			TString titleZ = Form("ResZL%d_R%d",layer,row);
			TString titleR = Form("ResRL%d_R%d",layer,row);
			TString titleW = Form("ResWL%d_R%d",layer,row);
			TString titleL = Form("ResLL%d_R%d",layer,row);
			histX[layer][row] = (TH1D*)file->Get(titleX);
			histY[layer][row] = (TH1D*)file->Get(titleY);
			histZ[layer][row] = (TH1D*)file->Get(titleZ);
			histR[layer][row] = (TH1D*)file->Get(titleR);
			histW[layer][row] = (TH1D*)file->Get(titleW);
			histL[layer][row] = (TH1D*)file->Get(titleL);
			if(tpc::Dead(layer,row) or histW[layer][row]->GetEffectiveEntries()<50){
				rx = 0;ry=0;rz = 0;rr=0;rw=0;rl=0;
				sx = 0;sy=0;sz = 0;sr=0;sw=0;sl=0;
			}
			else{
				peak = histX[layer][row]->GetBinCenter(histX[layer][row]->GetMaximumBin());
				width = histX[layer][row]->GetBinWidth(50);
				fGaus->SetRange(peak-4*width,peak+4*width);
				fGaus->SetParLimits(1,peak-2*width,peak+2*width);
				histX[layer][row] -> Fit("fGaus","QR");
				rx = fGaus->GetParameter(1);
				sx = fGaus->GetParameter(2);
			
				peak = histZ[layer][row]->GetBinCenter(histZ[layer][row]->GetMaximumBin());
				width = histZ[layer][row]->GetBinWidth(50);
				fGaus->SetRange(peak-3*width,peak+3*width);
				fGaus->SetParLimits(1,peak-2*width,peak+2*width);
				histZ[layer][row] -> Fit("fGaus","QR");
				rz = fGaus->GetParameter(1);
				sz = fGaus->GetParameter(2);
			
				peak = histR[layer][row]->GetBinCenter(histR[layer][row]->GetMaximumBin());
				width = histR[layer][row]->GetBinWidth(50);
				fGaus->SetRange(peak-4*width,peak+4*width);
				fGaus->SetParLimits(1,peak-2*width,peak+2*width);
				histR[layer][row] -> Fit("fGaus","QR");
				rr = fGaus->GetParameter(1);
				sr = fGaus->GetParameter(2);
				
				peak = histY[layer][row]->GetBinCenter(histY[layer][row]->GetMaximumBin());
				width = histY[layer][row]->GetBinWidth(50);
				fGaus->SetRange(peak-3*width,peak+3*width);
				fGaus->SetParLimits(1,peak-1*width,peak+1*width);
				histY[layer][row] -> Fit("fGaus","QR");
				ry = fGaus->GetParameter(1);
				sy = fGaus->GetParameter(2);
				
				peak = histW[layer][row]->GetBinCenter(histW[layer][row]->GetMaximumBin());
				width = histW[layer][row]->GetBinWidth(50);
				fGaus->SetRange(peak-4*width,peak+4*width);
				fGaus->SetParLimits(1,peak-2*width,peak+2*width);
				histW[layer][row] -> Fit("fGaus","QR");
				rw = fGaus->GetParameter(1);
				sw = fGaus->GetParameter(2);
				
				peak = histL[layer][row]->GetBinCenter(histL[layer][row]->GetMaximumBin());
				width = histL[layer][row]->GetBinWidth(50);
				fGaus->SetRange(peak-3*width,peak+3*width);
				fGaus->SetParLimits(1,peak-1*width,peak+1*width);
				histL[layer][row] -> Fit("fGaus","QR");
				rl = fGaus->GetParameter(1);
				sl = fGaus->GetParameter(2);
			}
			res_x.push_back(rx);
			res_y.push_back(ry);
			res_z.push_back(rz);
			res_r.push_back(rr);
			res_w.push_back(rw);
			res_l.push_back(rl);
			sig_x.push_back(sx);
			sig_y.push_back(sy);
			sig_z.push_back(sz);
			sig_r.push_back(sr);
			sig_w.push_back(sw);
			sig_l.push_back(sl);
		}
		bool Inner = true;
		if(layer > 9) Inner = false;
		vector<double> res_xCor = Smooth(res_x,0.5,Inner);
		vector<double> res_zCor = Smooth(res_z,0.5,Inner);
		vector<double> res_wCor = Smooth(res_w,0.5,Inner);
		AllResX.push_back(res_x);	
		AllResY.push_back(res_y);	
		AllResZ.push_back(res_z);	
		AllResW.push_back(res_w);	
		AllResL.push_back(res_l);	
		AllResR.push_back(res_r);	
		AllSigX.push_back(sig_x);	
		AllSigY.push_back(sig_y);	
		AllSigZ.push_back(sig_z);	
		AllSigR.push_back(sig_r);	
		AllSigW.push_back(sig_w);	
		AllSigL.push_back(sig_l);	
		AllResXCor.push_back(res_xCor);	
		AllResZCor.push_back(res_zCor);	
		AllResWCor.push_back(res_wCor);
	}
	cout<<"Writing..."<<endl;
	for(int layer =0;layer < 32;++layer){
		auto res_x = AllResX.at(layer);	
		auto res_y = AllResY.at(layer);	
		auto res_z = AllResZ.at(layer);	
		auto res_r = AllResR.at(layer);	
		auto res_w = AllResW.at(layer);	
		auto res_l = AllResL.at(layer);	
		auto sig_x = AllSigX.at(layer);	
		auto sig_y = AllSigY.at(layer);	
		auto sig_z = AllSigZ.at(layer);	
		auto sig_r = AllSigR.at(layer);	
		auto sig_w = AllSigW.at(layer);	
		auto sig_l = AllSigL.at(layer);	
		auto res_xCor = AllResXCor.at(layer);	
		auto res_zCor = AllResZCor.at(layer);	
		auto res_wCor = AllResWCor.at(layer);	
		int max_row = tpc::padParameter[layer][1];
		for(int row = 0;row < max_row; ++row){
			isDead = false;
			if(tpc::Dead(layer,row)) isDead = true;
			file2->cd();
			TString titleX = Form("ResXL%d_R%d",layer,row);
			TString titleY = Form("ResYL%d_R%d",layer,row);
			TString titleZ = Form("ResZL%d_R%d",layer,row);
			TString titleR = Form("ResRL%d_R%d",layer,row);
			TString titleW = Form("ResWL%d_R%d",layer,row);
			TString titleL = Form("ResLL%d_R%d",layer,row);
			histX[layer][row]->Write();
			histY[layer][row]->Write();
			histZ[layer][row]->Write();
			histR[layer][row]->Write();
			histW[layer][row]->Write();
			histL[layer][row]->Write();
			ent = histX[layer][row]->GetEffectiveEntries();l=layer;r=row;
			angle = tpc::getTheta(layer,row)*acos(-1)/180.;
			rx=res_x.at(row);
			ry=res_y.at(row);
			rz=res_z.at(row);
			rr=res_r.at(row);
			rw=res_w.at(row);
			rl=res_l.at(row);
			sx=sig_x.at(row);
			sy=sig_y.at(row);
			sz=sig_z.at(row);
			sr=sig_r.at(row);
			sw=sig_w.at(row);
			sl=sig_l.at(row);

			rxCor=res_xCor.at(row);
			rzCor=res_zCor.at(row);
			rwCor=res_wCor.at(row);
			CorX = - rw*cos(angle) + rl * sin(angle);
			CorZ =  rw*sin(angle) + rl*cos(angle);
			tree->Fill();
			double rxprev = -ParX[layer][row];
			double rzprev = -ParZ[layer][row];


			CorX +=rxprev;
			CorZ +=rzprev;
			bool correction = true;
			if(layer == 31 and row > 55) correction = true;
			if(layer == 30 and row > 65) correction = true;
			if(layer == 29 and row > 80) correction = true;
			if(!correction){
				CorX = 0;
				CorZ = 0;
			}
			//			f2<<Form("%d	%d	%d	%f	%f",layer,row,-250,-rxCor,-rzCor)<<endl;
			f2<<Form("%d	%d	%d	%g	%g",layer,row,-250,-CorX,-CorZ)<<endl;
		}
	}
	file2->Write();
}
void MakePadResidualHists(){
//	TFile* file = new TFile(Form("run0%d_DstTPCHelixTrackingHToF.root",runnum));
//	TFile* file = new TFile(Form("run0%d_DstTPCHelixTracking.root",runnum));
	TString filename = "./Beamthrough";
	TFile* file = new TFile(filename+".root");
	TTree* tree = (TTree*)file->Get("tpc");
	vector<double>*helix_cx = new vector<double>;
	vector<double>*helix_cy = new vector<double>;
	vector<double>*helix_r = new vector<double>;
	vector<double>*mom0 = new vector<double>;
	vector<double>*chisqr = new vector<double>;
	vector<double>*path = new vector<double>;
	vector<int>*helix_flag = new vector<int>;
	vector<int>*nhtrack = new vector<int>;
	vector<vector<double>>* hitpos_x = new vector<vector<double>>;
	vector<vector<double>>* hitpos_y = new vector<vector<double>>;
	vector<vector<double>>* hitpos_z = new vector<vector<double>>;
	vector<vector<double>>* residual_x = new vector<vector<double>>;
	vector<vector<double>>* residual_y = new vector<vector<double>>;
	vector<vector<double>>* residual_z = new vector<vector<double>>;
	vector<vector<double>>* theta_diff = new vector<vector<double>>;
	vector<vector<double>>* hitlayer = new vector<vector<double>>;
	vector<vector<double>>* track_cluster_row_center = new vector<vector<double>>;
	tree->SetBranchAddress("helix_cx",&helix_cx);
	tree->SetBranchAddress("helix_cy",&helix_cy);
	tree->SetBranchAddress("helix_r",&helix_r);
	tree->SetBranchAddress("helix_flag",&helix_flag);
	tree->SetBranchAddress("mom0",&mom0);
	tree->SetBranchAddress("chisqr",&chisqr);
	tree->SetBranchAddress("path",&path);
	tree->SetBranchAddress("nhtrack",&nhtrack);
	tree->SetBranchAddress("hitpos_x",&hitpos_x);
	tree->SetBranchAddress("hitpos_y",&hitpos_y);
	tree->SetBranchAddress("hitpos_z",&hitpos_z);
	tree->SetBranchAddress("residual_x",&residual_x);
	tree->SetBranchAddress("residual_y",&residual_y);
	tree->SetBranchAddress("residual_z",&residual_z);
	tree->SetBranchAddress("theta_diff",&theta_diff);
	tree->SetBranchAddress("hitlayer",&hitlayer);
	tree->SetBranchAddress("track_cluster_row_center",&track_cluster_row_center);
	

//	TFile* file2 = new TFile(Form("run0%d_ResidualHists.root",runnum),"create");
	TFile* file2 = new TFile(filename+"Hists.root","recreate");
	TH1D* histX[32][244];
	TH1D* histY[32][244];
	TH1D* histZ[32][244];
	TH1D* histR[32][244];
	TH1D* histL[32][244];
	TH1D* histW[32][244];
	TH1D* histLay[42];
	TH2D* histLayRes = new TH2D("ResR_Lay","ResR_Lay",42,-10,32,100,-5,5);
	for(int layer =0;layer < 32;++layer){
		int max_row = tpc::padParameter[layer][1]; 
		for(int row = 0;row < max_row; ++row){
			TString titleX = Form("ResXL%d_R%d",layer,row);
			TString titleY = Form("ResYL%d_R%d",layer,row);
			TString titleZ = Form("ResZL%d_R%d",layer,row);
			TString titleR = Form("ResRL%d_R%d",layer,row);
			TString titleW = Form("ResWL%d_R%d",layer,row);
			TString titleL = Form("ResLL%d_R%d",layer,row);
			histX[layer][row] = new TH1D(titleX,titleX,100,-15,15);	
			histY[layer][row] = new TH1D(titleY,titleY,100,-15,15);	
			histZ[layer][row] = new TH1D(titleZ,titleZ,100,-15,15);	
			histR[layer][row] = new TH1D(titleR,titleR,100,-15,15);	
			histW[layer][row] = new TH1D(titleW,titleW,100,-15,15);	
			histL[layer][row] = new TH1D(titleL,titleY,100,-15,15);	
		}
	}
	for(int layer = -10;layer < 32;++layer){
		TString title = Form("ResR_Layer%d",layer);
		histLay[layer+10] = new TH1D(title,title,100,-5,5);
	}
	int ent = tree->GetEntries();
	for(int iev = 0; iev < ent; ++iev){
//		if(iev%1==0) cout<<"Event "<<iev<<endl;
		tree->GetEntry(iev);
		int nt = residual_x->size();
		if(nt!=1) continue;
		for(int it=0;it<nt;++it){
			auto hf = helix_flag->at(it);
			if(hf != 100) continue;
			auto center_x = -helix_cx->at(it);
			auto center_z = helix_cy->at(it) - 143;
			auto r = helix_r->at(it);
			double target_dist = abs(hypot(center_x,center_z+143)-r);
//			if(target_dist>30) continue;
			auto hitx = hitpos_x->at(it);
			auto hity = hitpos_y->at(it);
			auto hitz = hitpos_z->at(it);
			auto resx = residual_x->at(it);
			auto resy = residual_y->at(it);
			auto resz = residual_z->at(it);
			auto hitlay = hitlayer->at(it);
			auto hitrow = track_cluster_row_center->at(it);
			int nh = resx.size();
			auto alpha = theta_diff->at(it);
//			if(nh < 15) continue;
			auto chi2 = chisqr -> at(it);
//			if(chi2>15) continue;
			auto pathl = path -> at(it);
//			if(pathl < 200 or pathl > 450) continue;
			auto p0 = mom0->at(it);
//			if(p0 - 0.10 < 0 or 1. -p0<0) continue;
			for(int ih=0;ih<nh;++ih){
//				if(abs(alpha.at(ih))>0.1) continue;
				double hitpy = hity.at(ih);
//				if(abs(hitpy)<50) continue;
				int layer = hitlay.at(ih);	
				int row = hitrow.at(ih);
				double angle = tpc::getTheta(layer,row)*acos(-1)/180.;
				double rx = resx.at(ih);
				double ry = resy.at(ih);
				double rz = resz.at(ih);
				double sin_ = sin(angle); 
				double cos_ = cos(angle);
				double rw = sin_*rz-cos_*rx;
				double rl = cos_*rz+sin_*rx;
				histX[layer][row]->Fill(rx);
				histY[layer][row]->Fill(ry);
				histZ[layer][row]->Fill(rz);
				histW[layer][row]->Fill(rw);
				histL[layer][row]->Fill(rl);
				double dx = center_x - hitx.at(ih); 
				double dz = center_z - hitz.at(ih); 
				double dr = hypot(dx,dz);
				double rr = dr - r;
				histR[layer][row]->Fill(rr);
				int ll;
				if(hitz.at(ih)<-143 and layer < 10) ll = -layer -1;
				else ll = layer;
				histLay[ll+10]->Fill(rr);
				histLayRes->Fill(ll,rr);
			}//ih
		}//it
		if(iev%10000==0) cout<<"Event "<<iev<<endl;
	}//iev
	file2->Write();
}



