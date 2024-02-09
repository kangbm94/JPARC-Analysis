void Sort(){
	TFile* file = new TFile("TPCInvM_cd15_CF1.root");
//	TFile* file = new TFile("old/TPCInv.root");
	TTree* tree = (TTree*)file->Get("tree");
	TFile* Sorted = new TFile("Sorted.root","recreate");
	TTree* SortedTree = new TTree("tree","tree");
	int runnum,evnum;
	bool FlgLd,FlgXi;
	double MM,InvMXi,mmpi0;
	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("evnum",&evnum);
	tree->SetBranchAddress("FlgLd",&FlgLd);
	tree->SetBranchAddress("FlgXi",&FlgXi);
	tree->SetBranchAddress("MM",&MM);
	tree->SetBranchAddress("mmpi0",&mmpi0);
	tree->SetBranchAddress("InvMXi",&InvMXi);
	SortedTree->Branch("runnum",&runnum);
	SortedTree->Branch("evnum",&evnum);
	int ent = tree->GetEntries();
	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
//		bool Acpt = abs(MM-1.315)<0.1 and mmpi0 > 0;// and abs(InvMXi-1.315)<0.1;
		bool Acpt =1;
//		if(FlgLd and FlgXi and Acpt){
			SortedTree->Fill();
			cout<<Form("(Run,ev) = (%d,%d)",runnum,evnum)<<endl;
//		}
	}
	cout<<SortedTree->GetEntries()<<endl;
	Sorted->Write();
}
