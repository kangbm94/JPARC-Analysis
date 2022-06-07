void test(){
	TFile* file = new TFile("dumb.root","RECREATE");
	TObjString s("comment1");
	file->WriteObject(&s,"test");
	file->Write();
}	
