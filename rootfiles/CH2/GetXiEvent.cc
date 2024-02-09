void GEtXiEvent(){
}
void GetXiEventKBM(){
//	TFile* f = new TFile("run05641_XiCand.root");
	TFile* f = new TFile("run05641_DstHSKKAna.root");
//	TFile* f = new TFile("run05641_DstE42.root");
	TTree* tree = (TTree*)f->Get("kk");
//	TTree* tree = (TTree*)f->Get("tpc");
	int evnum;
	bool Xiflag[5];
	bool KKflag[5];
/* vector<double>* MissMassCor=new vector<double>;
	vector<double>*pKurama=new vector<double>;
	vector<double>*qKurama=new vector<double>;
	vector<double>*chisqrKurama=new vector<double>;
	vector<double>*m2=new vector<double>;
	vector<int>* inside;
*/
	double MissMassCor[5],pKurama[5],qKurama[5],m2[5],chisqrKurama[5];
	tree->SetBranchAddress("evnum",&evnum);
//	tree->SetBranchAddress("MissMassCorr",&MissMassCor);
//	tree->SetBranchAddress("inside",inside);
	tree->SetBranchAddress("m2",m2);
	tree->SetBranchAddress("Xiflag",Xiflag);
	tree->SetBranchAddress("KKflag",KKflag);
	tree->SetBranchAddress("pKurama",pKurama);
	tree->SetBranchAddress("qKurama",qKurama);
	tree->SetBranchAddress("chisqrKurama",chisqrKurama);
//	double MissMassCord[5],pKuramad[5],qKuramad[5],m2d[5],chisqrKuramad[5];

	//	tree->SetBranchAddress("KKflag",&Xiflag);
//	tree->SetBranchAddress("Xiflag",&Xiflag);
//	tree->SetBranchAddress("MissMassCorr",MissMassCor);
//	TFile* f2 = new TFile("run05641_XiEvents_JWS.root","recreate");
	TFile* f2 = new TFile("run05641_XiEvents_KBM.root","recreate");
	TTree* tree2 = new TTree("tree","tree");
	tree2->Branch("evnum",&evnum);
	tree2->Branch("MissMassCorr",MissMassCor);
	for(int i=0;i<tree->GetEntries();++i){
		Xiflag[0]=false;
		tree->GetEntry(i);
//		cout<<i<<endl;
//		if(m2->size()==0) continue;
//		if(inside->size()==0) continue;
//		MissMassCord[0]=MissMassCor->at(0);
//		if(inside->at(0) and m2->at(0)>0.13 and m2->at(0)<0.4 and chisqrKurama->at(0)<200 and pKurama->at(0)<1.4 and qKurama->at(0)==1) Xiflag = true;	
		if(Xiflag[0]==true)tree2->Fill();
	}
	f2->Write();
	cout<<tree2->GetEntries()<<endl;
}
void GetXiEventJWS(){
	TFile* f = new TFile("run05641_DstE42.root");
	TTree* tree = (TTree*)f->Get("tpc");
	int evnum;
  vector<double>* MissMassCor=new vector<double>;
	vector<double>*pKurama=new vector<double>;
	vector<double>*qKurama=new vector<double>;
	vector<double>*chisqrKurama=new vector<double>;
	vector<double>*m2=new vector<double>;
	vector<int>* inside = new vector<int>;
	int nKK;
	tree->SetBranchAddress("evnum",&evnum);
	tree->SetBranchAddress("nKK",&nKK);
	tree->SetBranchAddress("MissMassCorr",&MissMassCor);
	tree->SetBranchAddress("inside",&inside);
	tree->SetBranchAddress("m2",&m2);
	tree->SetBranchAddress("pKurama",&pKurama);
	tree->SetBranchAddress("qKurama",&qKurama);
	tree->SetBranchAddress("chisqrKurama",&chisqrKurama);
	double MissMassCord[5],pKuramad[5],qKuramad[5],m2d[5],chisqrKuramad[5];
	bool Xiflag;
	//	tree->SetBranchAddress("KKflag",&Xiflag);
//	tree->SetBranchAddress("Xiflag",&Xiflag);
//	tree->SetBranchAddress("MissMassCorr",MissMassCor);
//	TFile* f2 = new TFile("run05641_XiEvents_JWS.root","recreate");
	TFile* f2 = new TFile("run05641_XiEvents_JWS.root","recreate");
	TTree* tree2 = new TTree("tree","tree");
	tree2->Branch("evnum",&evnum);
	tree2->Branch("MissMassCorr",MissMassCord);
	for(int i=0;i<tree->GetEntries();++i){
		Xiflag=false;
		tree->GetEntry(i);
//		cout<<i<<endl;
		if(nKK!= 1)continue;
		if(m2->size()==0) continue;
		if(inside->size()==0) continue;
//		MissMassCord[0]=MissMassCor->at(0);
		if(inside->at(0) and m2->at(0)>0.13 and m2->at(0)<0.4 and pKurama->at(0)<1.4 and qKurama->at(0)==1) Xiflag = true;	
		if(Xiflag==true)tree2->Fill();
	}
	f2->Write();
	cout<<tree2->GetEntries()<<endl;
}
