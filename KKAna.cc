#include "KKMethod.hh"
#include "KKLegacy.hh"
void KKAna(){
		KM.LoadFile("./rootfiles/CH2/AllKKAna.root");
//	KM.LoadFile("./rootfiles/CH2/DstKKAna05641.root");
	//		KM.LoadFile("../Other/E07_data/DstKKAna_CH2_Phase2.root");
	//	KM.LoadFile("~/WS_data/ch2target/run05666_KKAnaTest.root");
	KM.LoadKK();
}



void CountKPlus(){
	gStyle->SetOptFit(111);
	int nb=100;
	double m1= 0.35, m2 = 0.6;
	TChain* chain = KM.GetPublicChain();
	double KaonM2 = (KaonMass*KaonMass)/1000/1000;
	TH1D* hist = new TH1D("Kaons","Kaons",nb,m1,m2);
	chain->Draw("sqrt(m2)>>Kaons","pKurama<1.4&&inside==1&&qKurama>0");
	TF1* func = new TF1("GausWithPol","GausWithPol",0.3,0.7,5);
	TF1* func2 = new TF1("GausWithPol2","GausWithPol",0.3,0.7,5);
	double peak = hist->GetMaximum();
	double base = (hist->GetBinContent(1)+hist->GetBinContent(nb-1))/2;
	func->SetParLimits(0,peak/2,peak*2);
	func->SetParLimits(1,0.45,0.55);
	func->SetParLimits(2,0.01,0.05);
	func->SetRange(0.35,0.6);
	func->SetParLimits(3,base/2,base*2);
	hist->Fit("GausWithPol","R");
	double kpeak = func->GetParameter(0);
	double kmass = func->GetParameter(1);
	double ksig = func->GetParameter(2);
	double p1 = func->GetParameter(3);
	double p2 = func->GetParameter(4);
	double kparset[5]={kpeak,kmass,ksig,0,0};
	double bparset[5]={0,kmass,ksig,p1,p2};
	double nKaon = kpeak*sqrt(2*Pi())*ksig*nb/(m2-m1);
	func->SetParameters(kparset);
	func->SetLineColor(kBlue);
	func->Draw("same");

	func2->SetParameters(bparset);
	func2->SetLineColor(kBlack);
	func2->Draw("same");
	cout<<"Effective Entries: "<<hist->GetEffectiveEntries()<<endl;
	cout<<"Kaon Counts: "<<nKaon<<endl;
	cout<<"Background Counts: "<<(p1+p2*(m1/2+m2/2-kmass))*nb<<endl;

}

void GetPhi(){
	TChain* chain = KM.GetPublicChain();
	int ent = chain->GetEntries();
	cout<<ent<<endl;
	//	double p_cut=1.4;
	double KaonM2 = (KaonMass*KaonMass)/1000/1000;
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	TCanvas* c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(2,2);

	TFile* file = new TFile("KaonSorted.root","recreate");
	TTree* tree = new TTree("tree","tree");


	tree->Branch("nKaon",&nKaon,"nKaon/I");
	tree->Branch("BeamPx",BeamPx,"Beampx[nKaon]/D");
	tree->Branch("BeamPy",BeamPy,"Beampy[nKaon]/D");
	tree->Branch("BeamPz",BeamPz,"Beampz[nKaon]/D");
	tree->Branch("BeamP",BeamP,"Beampz[nKaon]/D");
	tree->Branch("KPx",KPx,"KPx[nKaon]/D");
	tree->Branch("KPy",KPy,"KPy[nKaon]/D");
	tree->Branch("KPz",KPz,"KPz[nKaon]/D");
	tree->Branch("KP",KP,"KP[nKaon]/D");
	tree->Branch("KM2",KM2,"KM2[nKaon]/D");
	tree->Branch("KQ",KQ,"KQ[nKaon]/I");

	ent=ent;
	TH2D* AllKaonPlot=new TH2D("Scat_AllKaons","Scat_AllKaons",100,-1,1,100,0,3);
	TH2D* TwoKaonPlot=new TH2D("Scat_TwoKaons","Scat_TwoKaons",100,-1,1,100,0,3);
	TH2D* KpTaggedPlot=new TH2D("KmScat_KpTagged","KmScat_KpTagged",100,-1,1,100,0,3);
	TH2D* KmTaggedPlot=new TH2D("KpScat_KmTagged","KpScat_KmTagged",100,-1,1,100,0,3);
	TH1D* PhiKPlot = new TH1D("Phi with Mk Fixed","Phi with Mk Fixed",100,0.95,1.35);
	TH1D* PhiPlot = new TH1D("Phi with Mk Not Fixed","Phi with Mk Not Fixed",100,0.95,1.35);
	bool flag[10];
	int nTotalTracks=0;
	int nTotalTracksChi2Cut=0;
	int nTotalTracksInside=0;
	int nDoubleKaons=0;
	int nTotalKaons=0;
	int nTotalKmKp=0;
	int nThreeOrMoreKaons=0;
	int nGoodKmKp=0;
	int TotalKp=0;
	int TotalKm=0;
	for(int i=0;i<ent;++i){
		Indicator(i,ent);
		chain->GetEntry(i);

		if(ntK18*ntKurama==0){
			continue;
		}
		nKaon=0;
		for(int ikk=0;ikk<25;ikk++){
			BeamPx[ikk]=-9999;	
			BeamPy[ikk]=-9999;	
			BeamPz[ikk]=-9999;	
			BeamP[ikk]=-9999;	
			
			KPx[ikk]=-9999;	
			KPy[ikk]=-9999;	
			KPz[ikk]=-9999;	
			KP[ikk]=-9999;
			KM2[ikk]=-9999;
			KQ[ikk]=-9999;
		} //		cout<<"Will Sort: "<<i<<endl; vector<Track> Tracks; Tracks.clear();
		vector<Track> Tracks;
		Tracks.clear();
		vector<Track> TracksChi2Cut;
		TracksChi2Cut.clear();
		vector<Track> TracksInside;
		TracksInside.clear();
		vector<Track> KaonTracks;
		KaonTracks.clear();
		int nkk=0;
		int nkkChi2Cut=0;
		int nkkInside=0;
		int nKaons=0;

		double K18Cut=20;
		double KuramaCut=200;
		for(int ikp=0;ikp<nKp;ikp++){
			for(int ikm=0;ikm<nKm;ikm++){
				int ikmkp=ikm+ikp*nKm;
				Track KKTrack(m2[ikp],pKurama[ikp],qKurama[ikp],ukp[ikmkp],vkp[ikmkp],inside[ikmkp],chisqrK18[ikm],chisqrKurama[ikp]);
				KKTrack.SetBeam(ukm[ikmkp],vkm[ikmkp],pK18[ikm]);
				Tracks.push_back(KKTrack);
			}
		}
		nkk=Tracks.size();
		nTotalTracks+=Tracks.size();
		bool inside_flag = false;
		for(int ikk=0;ikk<nkk;++ikk){
			if(Tracks[ikk].CutChiSqr(K18Cut,KuramaCut)){
				TracksChi2Cut.push_back(Tracks[ikk]);
				if(Tracks[ikk].IsInside()==1){
					inside_flag = true;
					if(Tracks[ikk].GetMomentum()<1.4) TracksInside.push_back(Tracks[ikk]);
				}
			}
			if(Tracks[ikk].CutChiSqr(K18Cut,KuramaCut)){
				if(Tracks[ikk].GetMomentum()<1.4&&inside_flag){
//					TracksInside.push_back(Tracks[ikk]);
				}
			}
		}
		nkkChi2Cut = TracksChi2Cut.size();
		nTotalTracksChi2Cut+=nkkChi2Cut;
		nkkInside = TracksInside.size();
		nTotalTracksInside+=nkkInside;

		for(int ikk=0;ikk<nkkInside;ikk++){
			if(abs(TracksInside[ikk].ParticleID())==KaonID){
				KaonTracks.push_back(TracksInside[ikk]);
			}
		}
		nKaons=KaonTracks.size();
		nTotalKaons+=nKaons;
		if(nKaons>2){
			nThreeOrMoreKaons++;
		}
		if(nKaons==2){
			nDoubleKaons++;
		}
		int ChargeSum = 0;
		for(int ikk=0;ikk<nKaons;++ikk){
			ChargeSum+=KaonTracks[ikk].GetCharge();
		}
		if(abs(ChargeSum)!=nKaons&&nKaons==2){
			nTotalKmKp++;
			double psum=KaonTracks[0].GetMomentum()+KaonTracks[1].GetMomentum();
			double psumx=KaonTracks[0].GetPx()+KaonTracks[1].GetPx();	
			double psumy=KaonTracks[0].GetPy()+KaonTracks[1].GetPy();	
			double psumz=KaonTracks[0].GetPz()+KaonTracks[1].GetPz();	
			double esumk = KaonTracks[0].GetKaonE()+KaonTracks[1].GetKaonE();	
			double esum = KaonTracks[0].GetE()+KaonTracks[1].GetE();	
			double psumMag = sqrt(psumx*psumx+psumy*psumy+psumz*psumz);
			double mphik = sqrt(esumk*esumk-psumMag);
			double mphi = sqrt(esum*esum-psumMag);
				nGoodKmKp++;	
				cout<<"PhiKMass : "<<mphik<<endl;
				cout<<"PhiMass : "<<mphi<<endl;
				PhiKPlot->Fill(mphik);
				PhiPlot->Fill(mphi);
		}

		for(int i=0;i<nKaons;i++){
			AllKaonPlot->Fill(KaonTracks[i].GetQM2(),KaonTracks[i].GetMomentum());
			if(nKaons==2){
				TwoKaonPlot->Fill(KaonTracks[i].GetQM2(),KaonTracks[i].GetMomentum());
			}
			if(KaonTracks[i].ParticleID()==KaonID&&KaonTracks[i].GetMomentum()<1.4){
				TotalKp++;
			}
			if(KaonTracks[i].ParticleID()==-KaonID&&KaonTracks[i].GetMomentum()<1.4){
				TotalKm++;
			}
		}
		if(nKaons==2){
			if(KaonTracks[0].ParticleID()>0){
				KpTaggedPlot->Fill(KaonTracks[1].GetQM2(),KaonTracks[1].GetMomentum());
			}
			else{
				KmTaggedPlot->Fill(KaonTracks[1].GetQM2(),KaonTracks[1].GetMomentum());
			}
			if(KaonTracks[1].ParticleID()>0){
				KpTaggedPlot->Fill(KaonTracks[0].GetQM2(),KaonTracks[0].GetMomentum());
			}
			else{
				KmTaggedPlot->Fill(KaonTracks[0].GetQM2(),KaonTracks[0].GetMomentum());
			}

		}
		nKaon=nKaons;
		for(int ikk=0;ikk<nKaon;ikk++){
			Track Kaon = KaonTracks[ikk];
			BeamPx[ikk]=Kaon.GetBeamPx();
			BeamPy[ikk]=Kaon.GetBeamPy();
			BeamPz[ikk]=Kaon.GetBeamPz();
			BeamP[ikk]=Kaon.GetBeamMomentum();
			
			KPx[ikk]=Kaon.GetPx();
			KPy[ikk]=Kaon.GetPy();
			KPz[ikk]=Kaon.GetPz();
			KP[ikk]=Kaon.GetMomentum();
			KM2[ikk]=Kaon.GetM2();
			KQ[ikk]=Kaon.GetCharge();
		}
		tree->Fill();
	}
	file->Write();
	cout<<"Total Tracks: "<<nTotalTracks<<endl;
	cout<<"Total Tracks From Vertex: "<<nTotalTracksInside<<endl;
	cout<<"Total Kaons: "<<nTotalKaons<<endl;
	cout<<"Total MultiKaons: "<<nThreeOrMoreKaons<<endl;
	cout<<"Total DoubleKaons: "<<nDoubleKaons<<endl;
	cout<<"Total KmKps "<<nTotalKmKp<<endl;
	cout<<"Total K+s "<<TotalKp<<endl;
	cout<<"Total K-s "<<TotalKm<<endl;
	cout<<"Total Good KmKps "<<nGoodKmKp<<endl;
	c1->cd(1);
	AllKaonPlot->Draw("colz");
	c1->cd(2);
	TwoKaonPlot->Draw("colz");
	c1->cd(3);
	KpTaggedPlot->Draw("colz");
	c1->cd(4);
	KmTaggedPlot->Draw("colz");
	c2->cd(1);
	PhiKPlot->Draw("colz");
	c2->cd(2);
	PhiPlot->Draw("colz");
}

void CountXi(double* par_,bool drawing){
	double par[6];
	TH1D* h = KM.XiMinusFit(par);
//	int nbin = h->GetXaxis()->GetNbinsX();
	int nbin = h->GetNbinsX();
	cout<<nbin<<endl;
	int bin_density = nbin/1;	
	TF1* XiSpectra = new TF1("XiSpectra","Gaussianf",1.2,1.4,3);
	TF1* XiBackground = new TF1("XiBackground","Gaussianf",1.2,1.4,3);
	XiSpectra->SetParameters(par[0],par[1],par[2]);
	XiSpectra->SetLineColor(kGreen);
	XiBackground->SetParameters(par[3],par[4],par[5]);
	XiBackground->SetLineColor(kBlue);
	GausWithBGf->SetParameters(par);
	GausWithBGf->SetLineColor(kRed);
	double XiCount = par[0]*par[2]*sqrt(2*Pi());
	cout<<"Number of Xis: "<<XiCount*bin_density<<endl;
	cout<<"Xi Mass: "<<1000*par[1]<<" MeV/c2"<<endl;
	cout<<"Xi Width: "<<1000*par[2]<<" MeV/c2"<<endl;
	for(int i=0;i<6;i++){
		par_[i]=par[i];
	}
	if(drawing){
		TCanvas* c1;
		c1= new TCanvas("XiCanv","XiCanv",1200,800);
		c1->cd();
		h->Draw();
		XiSpectra->Draw("same");
		XiBackground->Draw("same");
		GausWithBGf->Draw("same");
	}
}
void CountXi(){
	double par[6];bool draw = true;
	CountXi(par,draw);
}
void CountXiStar(double* par_,bool drawing){
	double par[6];
	TH1D* h = KM.XiStarFit(par);
//	int nbin = h->GetXaxis()->GetNbinsX();
	int nbin = h->GetNbinsX();
	cout<<nbin<<endl;
	int bin_density = nbin/1;	
	TF1* XiSpectra = new TF1("XiSpectra","Gaussianf",1.4,1.6,3);
	TF1* XiBackground = new TF1("XiBackground","Gaussianf",1.4,1.6,3);
	XiSpectra->SetParameters(par[0],par[1],par[2]);
	XiSpectra->SetLineColor(kGreen);
	XiBackground->SetParameters(par[3],par[4],par[5]);
	XiBackground->SetLineColor(kBlue);
	GausWithBGf->SetParameters(par);
	GausWithBGf->SetLineColor(kRed);
	double XiCount = par[0]*par[2]*sqrt(2*Pi());
	cout<<"Number of XiStars: "<<XiCount*bin_density<<endl;
	cout<<"Xi StarMass: "<<1000*par[1]<<" MeV/c2"<<endl;
	cout<<"Xi StarWidth: "<<1000*par[2]<<" MeV/c2"<<endl;
	for(int i=0;i<6;i++){
		par_[i]=par[i];
	}
	if(drawing){
		TCanvas* c1;
		c1= new TCanvas("XiCanv","XiCanv",1200,800);
		c1->cd();
		h->Draw();
		XiSpectra->Draw("same");
		XiBackground->Draw("same");
		GausWithBGf->Draw("same");
	}
}
void CountXiStar(){
	double par[6];bool draw = true;
	CountXiStar(par,draw);
}


void GetXiSpectra(){
	TChain* chain = KM.GetPublicChain();
	
	auto * h = KM.GetHistogram(6315);
	
	h->Draw();

	int ent = chain->GetEntries();
	cout<<ent<<endl;
	//	double p_cut=1.4;
	double KaonM2 = (KaonMass*KaonMass)/1000/1000;
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,2);
	TCanvas* c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(2,2);

	TFile* file = new TFile("KaonSorted.root","recreate");
	TTree* tree = new TTree("tree","tree");


	tree->Branch("nKaon",&nKaon,"nKaon/I");
	tree->Branch("BeamPx",BeamPx,"Beampx[nKaon]/D");
	tree->Branch("BeamPy",BeamPy,"Beampy[nKaon]/D");
	tree->Branch("BeamPz",BeamPz,"Beampz[nKaon]/D");
	tree->Branch("BeamP",BeamP,"Beampz[nKaon]/D");
	tree->Branch("KPx",KPx,"KPx[nKaon]/D");
	tree->Branch("KPy",KPy,"KPy[nKaon]/D");
	tree->Branch("KPz",KPz,"KPz[nKaon]/D");
	tree->Branch("KP",KP,"KP[nKaon]/D");
	tree->Branch("KM2",KM2,"KM2[nKaon]/D");
	tree->Branch("KQ",KQ,"KQ[nKaon]/I");

	ent=ent;
	TH2D* AllKaonPlot=new TH2D("Scat_AllKaons","Scat_AllKaons",100,-1,1,100,0,3);
	TH2D* TwoKaonPlot=new TH2D("Scat_TwoKaons","Scat_TwoKaons",100,-1,1,100,0,3);
	TH2D* KpTaggedPlot=new TH2D("KmScat_KpTagged","KmScat_KpTagged",100,-1,1,100,0,3);
	TH2D* KmTaggedPlot=new TH2D("KpScat_KmTagged","KpScat_KmTagged",100,-1,1,100,0,3);
	TH1D* PhiKPlot = new TH1D("Phi with Mk Fixed","Phi with Mk Fixed",100,0.95,1.35);
	TH1D* PhiPlot = new TH1D("Phi with Mk Not Fixed","Phi with Mk Not Fixed",100,0.95,1.35);
	bool flag[10];
	int nTotalTracks=0;
	int nTotalTracksChi2Cut=0;
	int nTotalTracksInside=0;
	int nDoubleKaons=0;
	int nTotalKaons=0;
	int nTotalKmKp=0;
	int nThreeOrMoreKaons=0;
	int nGoodKmKp=0;
	int TotalKp=0;
	int TotalKm=0;
	ent=ent/100;

	

	for(int i=0;i<ent;++i){
		Indicator(i,ent);
		chain->GetEntry(i);

		if(ntK18*ntKurama==0){
			continue;
		}
		nKaon=0;
		for(int ikk=0;ikk<25;ikk++){
			BeamPx[ikk]=-9999;	
			BeamPy[ikk]=-9999;	
			BeamPz[ikk]=-9999;	
			BeamP[ikk]=-9999;	
			
			KPx[ikk]=-9999;	
			KPy[ikk]=-9999;	
			KPz[ikk]=-9999;	
			KP[ikk]=-9999;
			KM2[ikk]=-9999;
			KQ[ikk]=-9999;
		} //		cout<<"Will Sort: "<<i<<endl; vector<Track> Tracks; Tracks.clear();
		vector<Track> Tracks;
		Tracks.clear();
		vector<Track> TracksChi2Cut;
		TracksChi2Cut.clear();
		vector<Track> TracksInside;
		TracksInside.clear();
		vector<Track> KaonTracks;
		KaonTracks.clear();
		int nkk=0;
		int nkkChi2Cut=0;
		int nkkInside=0;
		int nKaons=0;

		double K18Cut=20;
		double KuramaCut=200;
		for(int ikp=0;ikp<nKp;ikp++){
			for(int ikm=0;ikm<nKm;ikm++){
				int ikmkp=ikm+ikp*nKm;
				Track KKTrack(m2[ikp],pKurama[ikp],qKurama[ikp],ukp[ikmkp],vkp[ikmkp],inside[ikmkp],chisqrK18[ikm],chisqrKurama[ikp]);
				KKTrack.SetBeam(ukm[ikmkp],vkm[ikmkp],pK18[ikm]);
				Tracks.push_back(KKTrack);
			}
		}
		nkk=Tracks.size();
		nTotalTracks+=Tracks.size();
		bool inside_flag = false;
		for(int ikk=0;ikk<nkk;++ikk){
			if(Tracks[ikk].CutChiSqr(K18Cut,KuramaCut)){
				TracksChi2Cut.push_back(Tracks[ikk]);
				if(Tracks[ikk].IsInside()==1){
					inside_flag = true;
					if(Tracks[ikk].GetMomentum()<1.4) TracksInside.push_back(Tracks[ikk]);
				}
			}
			if(Tracks[ikk].CutChiSqr(K18Cut,KuramaCut)){
				if(Tracks[ikk].GetMomentum()<1.4&&inside_flag){
//					TracksInside.push_back(Tracks[ikk]);
				}
			}
		}
		nkkChi2Cut = TracksChi2Cut.size();
		nTotalTracksChi2Cut+=nkkChi2Cut;
		nkkInside = TracksInside.size();
		nTotalTracksInside+=nkkInside;

		for(int ikk=0;ikk<nkkInside;ikk++){
			if(abs(TracksInside[ikk].ParticleID())==KaonID){
				KaonTracks.push_back(TracksInside[ikk]);
			}
		}
		nKaons=KaonTracks.size();
		nTotalKaons+=nKaons;
		if(nKaons>2){
			nThreeOrMoreKaons++;
		}
		if(nKaons==2){
			nDoubleKaons++;
		}
		int ChargeSum = 0;
		for(int ikk=0;ikk<nKaons;++ikk){
			ChargeSum+=KaonTracks[ikk].GetCharge();
		}
		if(abs(ChargeSum)!=nKaons&&nKaons==2){
			nTotalKmKp++;
			double psum=KaonTracks[0].GetMomentum()+KaonTracks[1].GetMomentum();
			double psumx=KaonTracks[0].GetPx()+KaonTracks[1].GetPx();	
			double psumy=KaonTracks[0].GetPy()+KaonTracks[1].GetPy();	
			double psumz=KaonTracks[0].GetPz()+KaonTracks[1].GetPz();	
			double esumk = KaonTracks[0].GetKaonE()+KaonTracks[1].GetKaonE();	
			double esum = KaonTracks[0].GetE()+KaonTracks[1].GetE();	
			double psumMag = sqrt(psumx*psumx+psumy*psumy+psumz*psumz);
			double mphik = sqrt(esumk*esumk-psumMag);
			double mphi = sqrt(esum*esum-psumMag);
				nGoodKmKp++;	
				cout<<"PhiKMass : "<<mphik<<endl;
				cout<<"PhiMass : "<<mphi<<endl;
				PhiKPlot->Fill(mphik);
				PhiPlot->Fill(mphi);
		}

		for(int i=0;i<nKaons;i++){
			AllKaonPlot->Fill(KaonTracks[i].GetQM2(),KaonTracks[i].GetMomentum());
			if(nKaons==2){
				TwoKaonPlot->Fill(KaonTracks[i].GetQM2(),KaonTracks[i].GetMomentum());
			}
			if(KaonTracks[i].ParticleID()==KaonID&&KaonTracks[i].GetMomentum()<1.4){
				TotalKp++;
			}
			if(KaonTracks[i].ParticleID()==-KaonID&&KaonTracks[i].GetMomentum()<1.4){
				TotalKm++;
			}
		}
		if(nKaons==2){
			if(KaonTracks[0].ParticleID()>0){
				KpTaggedPlot->Fill(KaonTracks[1].GetQM2(),KaonTracks[1].GetMomentum());
			}
			else{
				KmTaggedPlot->Fill(KaonTracks[1].GetQM2(),KaonTracks[1].GetMomentum());
			}
			if(KaonTracks[1].ParticleID()>0){
				KpTaggedPlot->Fill(KaonTracks[0].GetQM2(),KaonTracks[0].GetMomentum());
			}
			else{
				KmTaggedPlot->Fill(KaonTracks[0].GetQM2(),KaonTracks[0].GetMomentum());
			}

		}
		nKaon=nKaons;
		for(int ikk=0;ikk<nKaon;ikk++){
			Track Kaon = KaonTracks[ikk];
			BeamPx[ikk]=Kaon.GetBeamPx();
			BeamPy[ikk]=Kaon.GetBeamPy();
			BeamPz[ikk]=Kaon.GetBeamPz();
			BeamP[ikk]=Kaon.GetBeamMomentum();
			
			KPx[ikk]=Kaon.GetPx();
			KPy[ikk]=Kaon.GetPy();
			KPz[ikk]=Kaon.GetPz();
			KP[ikk]=Kaon.GetMomentum();
			KM2[ikk]=Kaon.GetM2();
			KQ[ikk]=Kaon.GetCharge();
		}
		tree->Fill();
	}
	file->Write();
	cout<<"Total Tracks: "<<nTotalTracks<<endl;
	cout<<"Total Tracks From Vertex: "<<nTotalTracksInside<<endl;
	cout<<"Total Kaons: "<<nTotalKaons<<endl;
	cout<<"Total MultiKaons: "<<nThreeOrMoreKaons<<endl;
	cout<<"Total DoubleKaons: "<<nDoubleKaons<<endl;
	cout<<"Total KmKps "<<nTotalKmKp<<endl;
	cout<<"Total K+s "<<TotalKp<<endl;
	cout<<"Total K-s "<<TotalKm<<endl;
	cout<<"Total Good KmKps "<<nGoodKmKp<<endl;
	c1->cd(1);
	AllKaonPlot->Draw("colz");
	c1->cd(2);
	TwoKaonPlot->Draw("colz");
	c1->cd(3);
	KpTaggedPlot->Draw("colz");
	c1->cd(4);
	KmTaggedPlot->Draw("colz");
	c2->cd(1);
	PhiKPlot->Draw("colz");
	c2->cd(2);
	PhiPlot->Draw("colz");
}
