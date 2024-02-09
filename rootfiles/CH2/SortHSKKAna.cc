#include "HSKKAna.C"
void SortHSKKAna(){
//	TFile* file = new TFile("DstHSKKAna.root");
//	TTree* tree = (TTree*)file->Get("tree");
	TTree* tr = 0;	
	HSKKAna HK(tr);
	TFile* fileOut = new TFile("AllDstHSKKAnaXi.root","recreate");
	TTree* treeOut = new TTree("kk","Selected tree of HSKKAna");
	HK.AssignTree(treeOut);
	HK.AssignBranches();
	HK.ProcessTree();
	fileOut->Write();
}
