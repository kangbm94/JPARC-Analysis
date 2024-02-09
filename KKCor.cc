#include "KKMethod.hh"
#include "KKLegacy.hh"
#include "DetectorID.hh"
void KKCor(){
//		KM.LoadFile("./rootfiles/CH2/AllKKAna.root");
//		KM.LoadFile("./rootfiles/Production/AllKKAna.root");
//	KM.LoadFile("./rootfiles/CH2/DstKKAna05641.root");
	//		KM.LoadFile("../Other/E07_data/DstKKAna_CH2_Phase2.root");
	//	KM.LoadFile("~/WS_data/ch2target/run05666_KKAnaTest.root");
//	KM.LoadKK();
	gStyle->SetOptStat(10);
}

	
double dt = 0.025;


void GetXiSpectra(int runnum){
	TString infile,outfile;
	if(runnum==0){
		infile= "./rootfiles/CH2/DstHSKKAnaXi.root";
		outfile= "./SelectedEventsTest.root";
	}
	else{
		infile= "./rootfiles/CH2/run0"+to_string(runnum)+"_DstHSKKAna.root";
		outfile= "./KPTagged0"+to_string(runnum)+".root";
	}
	KM.LoadFile(infile);
	KM.LoadKK();
	TChain* chain = KM.GetPublicChain();
	KKEvent Event(chain);
	
	int ent = chain->GetEntries();
	int cnt =0;	
	TH1D* h = new TH1D("MissMass","MissMass",200,1,2);
	int XiEv,XiRun;
	double XiM2,XiP,XiU,XiV,XiTheta,XiThetaCM,KpMom;

	TFile* file = new TFile(outfile,"RECREATE"); TTree* tree = new TTree("tree","tree");
	tree->Branch("runnum",&XiRun,"runnum/I");	
	tree->Branch("evnum",&XiEv,"evnum/I");	
	tree->Branch("XiM2",&XiM2,"XiM2/D");	
	tree->Branch("XiMM",&XiM2,"XiMM/D");	
	tree->Branch("XiP",&XiP,"XiP/D");	
	tree->Branch("XiU",&XiU,"XiU/D");	
	tree->Branch("XiV",&XiV,"XiV/D");	
	tree->Branch("KpMom",&KpMom,"KpMom/D");	
	tree->Branch("XiTheta",&XiTheta,"XiTheta/D");	
	tree->Branch("XiThetaCM",&XiThetaCM,"XiThetaCM/D");	

	for(int i=0;i<ent;++i){
		Indicator(i,ent);
		Event.LoadEvent(i);
//		Event.CutTrig(22);//23 = DPS;
		Event.CutChiSqr();	
		Event.CutVertex();	
		Event.CutMomentum();	
		Event.CutCharge();	
		Event.CutM2();
		for(int ikk=0;ikk<Event.Getnkk();++ikk){
			if(Event.Getnkk()!=1) continue;
			XiRun=Event.GetRunNum();
			XiEv=Event.GetEvNum();
			XiM2=Event.GetMissMass(ikk);
			XiP=Event.GetMomentum(ikk);
			XiU=Event.GetUKP(ikk);
			XiV=Event.GetVKP(ikk);
			XiTheta=Event.GetTheta(ikk);
			XiThetaCM=Event.GetThetaCM(ikk);
			KpMom = Event.GetMomentum(ikk);
			cnt++;
			tree->Fill();
			h->Fill(XiM2);
			if(ikk){
				cout<<"Xi Multiplicity!"<<endl;
			}
		}
		Event.Clear();
	}
	file->Write();
	file->Close();
	cout<<cnt<<endl;
	h->Draw();
}
