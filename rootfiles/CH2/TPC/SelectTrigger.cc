vector<int>* trigflag = new vector<int>;
void SelectTrigger(){
	cout<<"SelectTrigger(TString filename, int trg)"<<endl;
}
void SelectTrigger(TString filename, int trg){
	TFile* file = new TFile(filename);
	TTree* tree = (TTree*)file->Get("tpc");
	tree->SetBranchAddress("trigflag",&trigflag);
	TFile* file2 = 
}
