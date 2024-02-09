double mK = 493.677/1000.;
double mP = 938.272/1000.;
double mXi = 1321.71/1000;
double KurCor = 1.;//;1.01453;
double gResolution[11];
auto fgaus = new TF1("fgaus","gaus",-0.2,0.2);
auto fpol4 = new TF1("fpol4","pol4",-0.2,0.2);
auto fpol5 = new TF1("fpol5","pol5",-0.4,0.15);

inline double XiMM(double p1,double u1,double v1,double p2,double u2,double v2){
	//p2 = p2*1.01453;
	double t1 = 1./sqrt(1+u1*u1+v1*v1);
	double p1z = p1*t1;
	double E1 = sqrt(mK*mK+p1*p1);
	double t2 = 1./sqrt(1+u2*u2+v2*v2);
	double p2z = p2*t2;
	p2 = p2;
	double E2 = sqrt(mK*mK+p2*p2);
	
	double E3 = E1+mP-E2;
	double p3 =sqrt(
		(u1*p1z-u2*p2z)*(u1*p1z-u2*p2z)+
		(v1*p1z-v2*p2z)*(v1*p1z-v2*p2z)+
		(p1z-p2z)*(p1z-p2z)
		);
	return sqrt(E3*E3-p3*p3);
}
vector<double>KmP; vector<double>KmU;
vector<double>KmV;
vector<double>KmX;
vector<double>KmY;
vector<double>KpP;
vector<double>KpU;
vector<double>KpV;
vector<double>KpX;
vector<double>KpY;
void fcnMM(int& npar, double* grad, double& fval, double* par,int flag){
	double du1 = par[0];
	double dv1 = par[1];
	double dp1 = par[2];
	double du2 = par[3];
	double dv2 = par[4];
	double dp2 = par[5];
	double lambda = par[6];
//	db2=1;
	int nh = KmP.size();
	fval = 0;
	double su1 = gResolution[0];
	double sv1 = gResolution[1];
	double sp1 = gResolution[2];
	double sx1 = gResolution[3];
	double sy1 = gResolution[4];
	double su2 = gResolution[5];
	double sv2 = gResolution[6];
	double sp2 = gResolution[7];
	double sx2 = gResolution[8];
	double sy2 = gResolution[9];
	double smm = gResolution[10];
	fval = 0;
	for(int i=0;i<nh;++i){
		double PKm = KmP.at(i)+dp1;
		double UKm = KmU.at(i)+du1;
		double VKm = KmV.at(i)+dv1;
		double PKp = KpP.at(i)+dp2;
		double UKp = KpU.at(i)+du2;
		double VKp = KpV.at(i)+dv2;
		double NKm = 1./sqrt(1+ UKm*UKm+ VKm*VKm);
		double NKp = 1./sqrt(1+ UKp*UKp+ VKp*VKp);
		double EKm = sqrt(mK*mK+PKm*PKm);	
		double EKp = sqrt(mK*mK+PKp*PKp);	
		double E = EKm+mP-EKp;
		double Px = PKm*UKm-PKp*UKp;
		double Py = PKm*VKm-PKp*VKp;
		double Pz = PKm*NKm-PKp*NKp;
	
		double P = sqrt(Px*Px+Py*Py+Pz*Pz);

		double mm = XiMM(PKm,UKm,VKm,PKp,UKp,VKp);
		
		double dPdu1 = (1./P)*( Px*PKm + Pz*NKm*NKm*NKm*UKm);
		double dPdv1 = (1./P)*( Py*PKm + Pz*NKm*NKm*NKm*VKm);
		double dPdp1 = (1./P)*( Px*UKm + Py*VKm + Pz*NKm);
		double dPdu2 = (1./P)*( -Px*PKp + -Pz*NKp*NKp*NKp*UKp);
		double dPdv2 = (1./P)*( -Py*PKp + -Pz*NKp*NKp*NKp*VKp);
		double dPdp2 = (1./P)*( -Px*UKm + -Py*VKm + -Pz*NKm);

		double dEdp1 = PKm/EKm;
		double dEdp2 = -PKp/EKp;

		double dXdu1 = du1/su1/su1 - lambda*(mm-mXi)*P*dPdu1/smm/smm/mm; 
		double dXdv1 = dv1/sv1/sv1 - lambda*(mm-mXi)*P*dPdv1/smm/smm/mm; 
		double dXdp1 = dp1/sp1/sp1 + lambda*(mm-mXi)*(E*PKm/EKm-P*dPdp1)/smm/smm/mm; 
		double dXdu2 = du2/su2/su2 - lambda*(mm-mXi)*P*dPdu2/smm/smm/mm; 
		double dXdv2 = dv2/sv2/sv2 - lambda*(mm-mXi)*P*dPdv2/smm/smm/mm; 
		double dXdp2 = dp2/sp2/sp2 + lambda*(mm-mXi)*(E*dEdp2-P*dPdp2)/smm/smm/mm; 
		double dXdL = (mm-mXi)*(mm-mXi)/smm/smm;

		fval += dXdu1*dXdu1 + dXdv1*dXdv1 + dXdp1*dXdp1
			+ dXdu2*dXdu2 + dXdv2*dXdv2 + dXdp2*dXdp2 + dXdL*dXdL;
	}
	par[7]=fval;
	fval = fval / nh;
}
void fcnVtx(int& npar, double* grad, double& fval, double* par,int flag){
	double xkm_cor = par[0];
	double ykm_cor = par[1];
	double ukm_cor = par[2];
	double vkm_cor = par[3];
	double xkp_cor = par[4];
	double ykp_cor = par[5];
	double ukp_cor = par[6];
	double vkp_cor = par[7];
	int nh = KmX.size();
	for(int i=0;i<nh;++i){
		double xkm = KmX.at(i)+xkm_cor;
		double ykm = KmY.at(i)+ykm_cor;
		double ukm = KmU.at(i)+ukm_cor;
		double vkm = KmV.at(i)+vkm_cor;
		double xkp = KpX.at(i)+xkp_cor;
		double ykp = KpY.at(i)+ykp_cor;
		double ukp = KpU.at(i)+ukp_cor;
		double vkp = KpV.at(i)+vkp_cor;
		double z = (xkm-xkp)*(ukp-ukm)+(ykm-ykp)*(vkp-vkm)/
		((ukp-ukm)*(ukp-ukm)+(vkp-vkm)*(vkp-vkm));
	

	}



}
void GetKMKP(){
	cout<<"CorrectKMKP"<<endl;
	cout<<"CorrectPCalc"<<endl;
	cout<<"CorrectVertZ"<<endl;
	cout<<"GetKMKPFile"<<endl;
	double vtxSig = 6.5855;
	double vtySig = 6.5855;
	double vtzSig = 118;
	double K18Sig = 0.0021;
	double K18USig = 0.00304;
	double K18VSig = 0.00268;
	double KuramaSig = 0.042;
	double KuramaUSig = 0.00304;
	double KuramaVSig = 0.00269;
	gResolution[0]=K18USig;
	gResolution[1]=K18VSig;
	gResolution[2]=K18Sig;
	gResolution[5]=KuramaUSig;
	gResolution[6]=KuramaVSig;
	gResolution[7]=KuramaSig;
	gResolution[10]=0.05;

}
void CorrectPCalc(){
//	TFile* file = new TFile("NoCorrection/KMKP.root");
	TFile* file = new TFile("KMKP.root");
	TTree* tree = (TTree*)file->Get("tree"); 
	double KMP,KMU,KMV,KPP,KPU,KPV,KpCalc;
	tree->SetBranchAddress("KMP",&KMP);
	tree->SetBranchAddress("KMU",&KMU);
	tree->SetBranchAddress("KMV",&KMV);
	tree->SetBranchAddress("KPU",&KPU);
	tree->SetBranchAddress("KPV",&KPV);
	tree->SetBranchAddress("KPP",&KPP);
	tree->SetBranchAddress("KpCalc",&KpCalc);
	int ent = tree->GetEntries();
	int nbX = 500;
	TH2D* PUM = new TH2D("PUM","PUM",nbX,-0.15,0.05,200,-1,1);
	TH2D* PUP = new TH2D("PUP","PUP",nbX,-0.3,0.05,100,-0.5,0.5);
	TH2D* PVP = new TH2D("PVP","PVP",nbX,-0.15,0.15,100,-0.5,0.5);
	int nslice = 20;
	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
		PUM->Fill(KMU,KPP-KpCalc);
		PUP->Fill(KPU,KPP-KpCalc);
		PVP->Fill(KPV,KPP-KpCalc);
	}
	TH1D* PUM_h[100];
	TH1D* PUP_h[100];
	TH1D* PVP_h[100];
	TCanvas* c1 = new TCanvas("c1","c1",800,800);
	c1->Divide(5,4);
	TCanvas* c2 = new TCanvas("c2","c2",800,800);
	c2->Divide(5,4);
	vector<double>UPC,VPC;
	TGraphErrors* PUPG = new TGraphErrors();
	TGraphErrors* PVPG = new TGraphErrors();
	for(int i = 0; i<nslice; ++i){
		int nbin = nbX / nslice;
		int b1 = i*nbin+1;
		int b2 = (i+1)*nbin;
		double up1 = (PUP->GetXaxis())->GetBinCenter(b1);
		double up2 = (PUP->GetXaxis())->GetBinCenter(b2);
		double up = (up1+up2)/2;
		double upw = (PUP->GetXaxis())->GetBinWidth(b2);
		double vp1 = (PVP->GetXaxis())->GetBinCenter(b1);
		double vp2 = (PVP->GetXaxis())->GetBinCenter(b2);
		double vp = (vp1+vp2)/2;
		double vpw = (PVP->GetXaxis())->GetBinWidth(b2);
		c1->cd(i+1);
		PUP_h[i] = PUP->ProjectionY(Form("PUP%d",i),nbin*i,nbin*(i+1));
		PUP_h[i]->Draw();
		int peakbin = PUP_h[i]->GetMaximumBin();
		double upc = PUP_h[i]->GetBinCenter(peakbin);
		fgaus->SetRange(upc-0.03,upc+0.03);
		PUP_h[i]->Fit("fgaus","R");	
		upc = fgaus->GetParameter(1);
		double ups = fgaus->GetParameter(2);
		UPC.push_back(upc);
		if((i != 0 and i < nslice-2) or 1){
			PUPG->AddPoint(up,upc);
			PUPG->SetPointError(PUPG->GetN()-1,upw/sqrt(12),ups);
		}
		c2->cd(i+1);
		PVP_h[i] = PVP->ProjectionY(Form("PVP%d",i),nbin*i,nbin*(i+1));
		PVP_h[i]->Draw();
		peakbin = PVP_h[i]->GetMaximumBin();
		double vpc = PVP_h[i]->GetBinCenter(peakbin);
		fgaus->SetRange(vpc-0.03,vpc+0.03);
		PVP_h[i]->Fit("fgaus","R");	
		vpc = fgaus->GetParameter(1);
		double vps = fgaus->GetParameter(2);
		PVPG->AddPoint(vp,vpc);
		PVPG->SetPointError(PVPG->GetN()-1,vpw/sqrt(12),vps);
		VPC.push_back(vpc);
	}
	TCanvas* c3 = new TCanvas("c3","c3",50,50,800,800);
	c3->Divide(2);
	c3->cd(1);
	PUP->Draw("colz");
	PUPG->Draw("same");
	PUPG->Fit("fpol5");
	cout<<"Correcting dP vs U"<<endl;
	for(int i=0;i<6;++i){
		cout<<fpol5->GetParameter(i)<<endl;
	}
	c3->cd(2);
	PVP->Draw("colz");
	PVPG->Draw("same");
	PVPG->Fit("fpol4");
	cout<<"Correcting dP vs V"<<endl;
	for(int i=0;i<5;++i){
		cout<<fpol4->GetParameter(i)<<endl;
	}
	c3->cd(2);

}
void CorrectVertZ(){
//	TFile* file = new TFile("KMKP.root");
	TFile* file = new TFile("NoVertCor/KMKP.root");
//	TFile* file = new TFile("BeforeVCor/KMKP.root");
	TTree* tree = (TTree*)file->Get("tree"); 
	double KMP,KPP;
	double KMU,KMV,KPU,KPV,vtz;
	int ent = tree->GetEntries();
	int nbX = 100;
	TH2D* VUP = new TH2D("VUP","VUP",nbX,-0.40,0.05,50,-140,140);
	TH2D* VVP = new TH2D("VVP","VVP",nbX,-0.15,0.15,50,-140,140);
	tree->SetBranchAddress("KPU",&KPU);
	tree->SetBranchAddress("KPV",&KPV);
	tree->SetBranchAddress("vtz",&vtz);
	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
		VUP->Fill(KPU,vtz);
		VVP->Fill(KPV,vtz);
	}
	TH1D* VUP_h[100];
	TH1D* VVP_h[100];
	TCanvas* c1 = new TCanvas("c1","c1",800,800);
	c1->Divide(5,4);
	TCanvas* c2 = new TCanvas("c2","c2",800,800);
	c2->Divide(5,4);
	int nslice = 20;	
	auto fgaus = new TF1("fgaus","gaus",-0.2,0.2);
	TGraphErrors* VUPG = new TGraphErrors();
	TGraphErrors* VVPG = new TGraphErrors();
	for(int i = 0; i<nslice; ++i){
		int nbin = nbX / nslice;
		int b1 = i*nbin+1;
		int b2 = (i+1)*nbin;
		double up1 = (VUP->GetXaxis())->GetBinCenter(b1);
		double up2 = (VUP->GetXaxis())->GetBinCenter(b2);
		double up = (up1+up2)/2;
		double upw = (VUP->GetXaxis())->GetBinWidth(b2);
		double vp1 = (VVP->GetXaxis())->GetBinCenter(b1);
		double vp2 = (VVP->GetXaxis())->GetBinCenter(b2);
		double vp = (vp1+vp2)/2;
		double vpw = (VVP->GetXaxis())->GetBinWidth(b2);
		c1->cd(i+1);
		VUP_h[i] = VUP->ProjectionY(Form("VUP%d",i),nbin*i,nbin*(i+1));
		VUP_h[i]->Draw();
		int peakbin = VUP_h[i]->GetMaximumBin();
		double upc = VUP_h[i]->GetBinCenter(peakbin);
		fgaus->SetRange(upc-30,upc+30);
		fgaus->SetParLimits(1,upc-10,upc+10);
		VUP_h[i]->Fit("fgaus","R");	
		upc = fgaus->GetParameter(1);
		double ups = fgaus->GetParameter(2);
		VUPG->AddPoint(up,upc);
		VUPG->SetPointError(VUPG->GetN()-1,upw/sqrt(12),ups);
		c2->cd(i+1);
		VVP_h[i] = VVP->ProjectionY(Form("VVP%d",i),nbin*i,nbin*(i+1));
		VVP_h[i]->Draw();
		peakbin = VVP_h[i]->GetMaximumBin();
		double vpc = VVP_h[i]->GetBinCenter(peakbin);
		fgaus->SetRange(vpc-30,vpc+30);
		fgaus->SetParLimits(1,vpc-10,vpc+10);
		VVP_h[i]->Fit("fgaus","R");	
		vpc = fgaus->GetParameter(1);
		double vps = fgaus->GetParameter(2);
		VVPG->AddPoint(vp,vpc);
		VVPG->SetPointError(VVPG->GetN()-1,vpw/sqrt(12),vps);
	
	}
	TCanvas* c3 = new TCanvas("c3","c3",50,50,800,800);
	c3->Divide(2);
	c3->cd(1);
	VUP->Draw("colz");
	VUPG->Draw("same");
	VUPG->Fit("fpol5");
	double p0,p1,p2,p3,p4,p5;
	p0 = fpol5->GetParameter(0);
	p1 = fpol5->GetParameter(1);
	p2 = fpol5->GetParameter(2);
	p3 = fpol5->GetParameter(3);
	p4 = fpol5->GetParameter(4);
	p5 = fpol5->GetParameter(5);
	cout<<Form("Param : (%f,%f,%f,%f,%f,%f)",p0,p1,p2,p3,p4,p5)<<endl;
	c3->cd(2);
	VVP->Draw("colz");
	VVPG->Draw("same");
	VVPG->Fit("fpol4");
	p0 = fpol4->GetParameter(0);
	p1 = fpol4->GetParameter(1);
	p2 = fpol4->GetParameter(2);
	p3 = fpol4->GetParameter(3);
	p4 = fpol4->GetParameter(4);
	cout<<Form("Param : (%f,%f,%f,%f,%f)",p0,p1,p2,p3,p4)<<endl;
}
void CorrectKMKP(){
	TFile* file = new TFile("KMKP.root");
	TTree* tree = (TTree*)file->Get("tree"); 
	double KMP,KPP;
	double KMU,KMV,KPU,KPV;
	tree->SetBranchAddress("KMP",&KMP);
	tree->SetBranchAddress("KMU",&KMU);
	tree->SetBranchAddress("KMV",&KMV);
	tree->SetBranchAddress("KPP",&KPP);
	tree->SetBranchAddress("KPU",&KPU);
	tree->SetBranchAddress("KPV",&KPV);
	int ent = tree->GetEntries();
	TH2D* MMU = new TH2D("MMU","U : mXi",50,-0.5,0.1,50,mXi-0.15,mXi+0.15);
	TH2D* MMV = new TH2D("MMV","V : mXi",50,-0.3,0.2,50,mXi-0.15,mXi+0.15);
	TH2D* MMP = new TH2D("MMP","P : mXi",50,1.75,1.9,50,mXi-0.15,mXi+0.15);
	TH2D* MMP2 = new TH2D("MMP2","P2 : mXi",50,1,1.5,50,mXi-0.15,mXi+0.15);
	TH2D* MMUCor = new TH2D("MMUCor","UCor : mXi",50,-0.5,0.1,50,mXi-0.1,mXi+0.1);
	TH2D* MMVCor = new TH2D("MMVCor","VCor : mXi",50,-0.3,0.2,50,mXi-0.1,mXi+0.1);
	TH2D* MMPCor = new TH2D("MMP","P : mXi",50,1.75,1.9,50,mXi-0.1,mXi+0.1);
	TH2D* MMP2Cor = new TH2D("MMP2CorOff","P2Cor : mXi",50,1,1.5,50,mXi-0.1,mXi+0.1);
	TH2D* MMUOff = new TH2D("MMUOff","U : mXi",50,-0.5,0.1,50,mXi-0.15,mXi+0.15);
	TH2D* MMVOff = new TH2D("MMVOff","V : mXi",50,-0.3,0.2,50,mXi-0.15,mXi+0.15);
	TH2D* MMPOff = new TH2D("MMPOff","P : mXi",50,1.75,1.9,50,mXi-0.15,mXi+0.15);
	TH2D* MMP2Off = new TH2D("MMP2Off","P2 : mXi",50,1,1.5,50,mXi-0.15,mXi+0.15);
	TH1D* HU1 = new TH1D("U1Cor","U1Cor",1000,-0.001,0.001);
	TH1D* HV1 = new TH1D("V1Cor","V1Cor",1000,-0.001,0.001);
	TH1D* HP1 = new TH1D("P1Cor","P1Cor",1000,-20,20);
	TH1D* HU2 = new TH1D("U2Cor","U2Cor",1000,-0.001,0.001);
	TH1D* HV2 = new TH1D("V2Cor","V2Cor",1000,-0.001,0.001);
	TH1D* HP2 = new TH1D("P2Cor","P2Cor",1000,-20,20);


	TMinuit* min = new TMinuit(6);
	min->SetFCN(fcnMM);
	double u1cor,u1cer;
	double v1cor,v1cer;
	double p1cor,p1cer;
	double u2cor,u2cer;
	double v2cor,v2cer;
	double p2cor,p2cer;
	double arglist[30];
	int ierflg = 0;
	arglist[0]=1;
	double step = 0.0001;
	double pstep = 0.1;
	arglist[0]=1000;//maxcalls
	arglist[1]=0.00001;//tolerance
	double du1,dv1,dp1,du2,dv2,dp2;
	double Xsum = 0;
	double X,Xcer;
//	gResolution[0]+= 0.003;
	for(int i=0;i<ent; ++i){
		tree->GetEntry(i);
//		KPU+=0.1;	
//		KMU= KMU+0.005;
		KmP.push_back(KMP);
		KmU.push_back(KMU);
		KmV.push_back(KMV);
		KpP.push_back(KPP);
		KpU.push_back(KPU);
		KpV.push_back(KPV);
		double mm = XiMM(KMP,KMU,KMV,KPP*KurCor,KPU,KPV);
		MMU->Fill(KPU,mm);
		MMV->Fill(KPV,mm);
		MMP->Fill(KMP,mm);
		MMP2->Fill(KPP,mm);
		min->mnparm(0,"u1cor",0.,step,-0.01,0.01,ierflg);
		min->mnparm(1,"v1cor",0.,step,-0.01,0.01,ierflg);
		min->mnparm(2,"p1cor",0.,pstep,-20,20,ierflg);
		min->mnparm(3,"u2cor",0.,step,-0.01,0.01,ierflg);
		min->mnparm(4,"v2cor",0,step,-0.01,0.01,ierflg);
		min->mnparm(5,"p2cor",0,pstep,-20,20,ierflg);
		min->mnexcm("MINIMIZE",arglist,6,ierflg);
		min->GetParameter(0,u1cor,u1cer);
		min->GetParameter(1,v1cor,v1cer);
		min->GetParameter(2,p1cor,p1cer);
		min->GetParameter(3,u2cor,u2cer);
		min->GetParameter(4,v2cor,v2cer);
		min->GetParameter(5,p2cor,p2cer);
		min->GetParameter(7,X,Xcer);
		Xsum += 1;
		KmP.clear();
		KmU.clear();
		KmV.clear();
		KpP.clear();
		KpU.clear();
		KpV.clear();
		du1+=u1cor;	
		dv1+=v1cor;	
		dp1+=p1cor;	
		du2+=u2cor;	
		dv2+=v2cor;	
		dp2+=p2cor;	
//		HU1->Fill(u1cor);	
//		cout<<u1cor<<endl;
//		HV1->Fill(v1cor);	
//		HP1->Fill(p1cor);	
//		HU2->Fill(u2cor);	
//		HV2->Fill(v2cor);	
//		HP2->Fill(p2cor);	
		mm = XiMM(KMP+p1cor,KMU+u1cor,KMV+v1cor,KPP+p2cor,KPU+u2cor,KPV+v2cor);
		MMUCor->Fill(KPU,mm);
		MMVCor->Fill(KPV,mm);
		MMPCor->Fill(KMP,mm);
		MMP2Cor->Fill(KPP,mm);

	}
	du1 /= Xsum;	
	dv1 /= Xsum;	
	dp1 /= Xsum;	
	du2 /= Xsum;	
	dv2 /= Xsum;	
	dp2 /= Xsum;	
	cout<<Form("params : (%f,%f,%f,%f,%f,%f)",du1,dv1,dp1,du2,dv2,dp2)<<endl;
	

	for(int i=0;i<ent; ++i){
		tree->GetEntry(i);
		double mm = XiMM(KMP+dp1,KMU+du1,KMV+dv1,KPP+dp2,KPU+du2,KPV+dv2);
		MMUOff->Fill(KPU,mm);
		MMVOff->Fill(KPV,mm);
		MMPOff->Fill(KMP,mm);
		MMP2Off->Fill(KPP,mm);
//		cout<<mm<<endl;
	}
	TCanvas*c1 = new TCanvas("c1","c1",900,900);
	c1->Divide(4,3);
	c1->cd(1);
	MMU->Draw("colz");
	c1->cd(2);
	MMV->Draw("colz");
	c1->cd(3);
	MMP->Draw("colz");
	c1->cd(4);
	MMP2->Draw("colz");
	c1->cd(5);
	MMUCor->Draw("colz");
	c1->cd(6);
	MMVCor->Draw("colz");
	c1->cd(7);
	MMPCor->Draw("colz");
	c1->cd(8);
	MMP2Cor->Draw("colz");
	c1->cd(9);
	MMUOff->Draw("colz");
	c1->cd(10);
	MMVOff->Draw("colz");
	c1->cd(11);
	MMPOff->Draw("colz");
	c1->cd(12);
	MMP2Off->Draw("colz");
#if 0
	TCanvas*c2 = new TCanvas("c2","c2",900,900);
	c2->Divide(3,2);
	c2->cd(1);
	HU1->Draw();
	c2->cd(2);
	HV1->Draw();
	c2->cd(3);
	HP1->Draw();
	c2->cd(4);
	HU2->Draw();
	c2->cd(5);
	HV2->Draw();
	c2->cd(6);
	HP2->Draw();
#endif
}
void GetKMKPFile(){
	TFile* file = new TFile("AllDstHSKKAnaXi.root");
	TTree* tree = (TTree*)file->Get("kk");
	TFile* file2 = new TFile("KMKP.root","recreate");
	TTree* tree2 = new TTree("tree","tree");
	int nKK;
	double KMPX[5],KMPY[5],KMPZ[5];
	double KPPX[5],KPPY[5],KPPZ[5];
	double vtx[5],vty[5],vtz[5];
	double xtgtKurama[5],ytgtKurama[5],ztgtKurama[5];
	double xtgtHS[5],ytgtHS[5],ztgtHS[5];
	double ukp[5],vkp[5],pCalc[5],pCorr[5];
	bool Xiflag[5];
	double KMP,KPP;
	double KMU,KMV,KPU,KPV;
	double KpCalc,mm; 
	tree->SetBranchAddress("xtgtKurama",xtgtKurama);
	tree->SetBranchAddress("ytgtKurama",ytgtKurama);
	tree->SetBranchAddress("ztgtKurama",ztgtKurama);
	tree->SetBranchAddress("xtgtHS",xtgtHS);
	tree->SetBranchAddress("ytgtHS",ytgtHS);
	tree->SetBranchAddress("ztgtHS",ztgtHS);
	tree->SetBranchAddress("nKK",&nKK);

	tree->SetBranchAddress("ukp",ukp);
	tree->SetBranchAddress("vkp",vkp);
	tree->SetBranchAddress("vtx",vtx);
	tree->SetBranchAddress("vty",vty);
	tree->SetBranchAddress("vtz",vtz);
	tree->SetBranchAddress("pCalc",pCalc);
	tree->SetBranchAddress("pCorr",pCorr);

	tree->SetBranchAddress("KMPX",KMPX);
	tree->SetBranchAddress("KMPY",KMPY);
	tree->SetBranchAddress("KMPZ",KMPZ);
	tree->SetBranchAddress("KPPX",KPPX);
	tree->SetBranchAddress("KPPY",KPPY);
	tree->SetBranchAddress("KPPZ",KPPZ);
	tree->SetBranchAddress("Xiflag[nKK]",Xiflag);
	int ent = tree->GetEntries();
	tree2->Branch("nKK",&nKK);
	tree2->Branch("KMP",&KMP);
	tree2->Branch("KMU",&KMU);
	tree2->Branch("KMV",&KMV);
	tree2->Branch("KPP",&KPP);
	tree2->Branch("KPU",&KPU);
	tree2->Branch("KPV",&KPV);
	tree2->Branch("KpCalc",&KpCalc);
	tree2->Branch("mm",&mm);
	tree2->Branch("vtx",&vtx[0]);
	tree2->Branch("vty",&vty[0]);
	tree2->Branch("vtz",&vtz[0]);
	tree2->Branch("xtgtKurama",&xtgtKurama[0]);
	tree2->Branch("ytgtKurama",&ytgtKurama[0]);
	tree2->Branch("ztgtKurama",&ztgtKurama[0]);
	tree2->Branch("xtgtHS",&xtgtHS[0]);
	tree2->Branch("ytgtHS",&ytgtHS[0]);
	tree2->Branch("ztgtHS",&ztgtHS[0]);
	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
//		if(nKK!=1)continue;
//		if(!Xiflag[0]) continue;
		KMP = sqrt(KMPX[0]*KMPX[0]+KMPY[0]*KMPY[0]+KMPZ[0]*KMPZ[0]);
		KPP = sqrt(KPPX[0]*KPPX[0]+KPPY[0]*KPPY[0]+KPPZ[0]*KPPZ[0]);
		KMU = KMPX[0]/KMPZ[0];
		KMV = KMPY[0]/KMPZ[0];
		KPU = KPPX[0]/KPPZ[0];
		KPV = KPPY[0]/KPPZ[0];
		KpCalc = pCalc[0];
		if(KPU != ukp[0]){
			cout<<"Warning"<<endl;
		}
		mm = XiMM(KMP,KMU,KMV,KPP,KPU,KPV);
//		if(abs(mm-mXi)>0.02)	continue;
		tree2->Fill();
	}
	file2->Write();
}
