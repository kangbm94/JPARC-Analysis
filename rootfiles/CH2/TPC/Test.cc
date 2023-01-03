int ntTpc,nclTpc,runnum,evnum,process_time;
vector<int>* pid = new vector<int>;
vector<int>* charge = new vector<int>;
vector<int>* isBeam = new vector<int>;
vector<double>* mom0 = new vector<double>;
vector<double>* mom_vtx = new vector<double>;
vector<double>* mom_vty = new vector<double>;
vector<double>* mom_vtz = new vector<double>;
vector<double>* helix_cx = new vector<double>;
vector<double>* helix_cy = new vector<double>;
vector<double>* helix_z0 = new vector<double>;
vector<double>* helix_r = new vector<double>;
vector<double>* helix_dz = new vector<double>;
vector<double>* closeDist = new vector<double>;
vector<double>* chisqr = new vector<double>;
vector<double>* vtx = new vector<double>;
vector<double>* vty = new vector<double>;
vector<double>* vtz = new vector<double>;
vector<double>* cluster_de = new vector<double>;
vector<double>* cluster_x = new vector<double>;
vector<double>* cluster_y = new vector<double>;
vector<double>* cluster_z = new vector<double>;
vector<int>* cluster_size = new vector<int>;
vector<int>* hough_flag = new vector<int>;
vector<double>* tcluster_x = new vector<double>;
vector<double>* tcluster_y = new vector<double>;
vector<double>* tcluster_z = new vector<double>;
vector<int>* tcluster_size = new vector<int>;
vector<int>* combi_id = new vector<int>;
vector<vector<double>>* track_cluster_size = new vector<vector<double>>;
vector<vector<double>>* track_cluster_de = new vector<vector<double>>;
vector<vector<double>>* track_cluster_x_center = new vector<vector<double>>;
vector<vector<double>>* track_cluster_y_center = new vector<vector<double>>;
vector<vector<double>>* track_cluster_z_center = new vector<vector<double>>;
int helixid=0;
double PI = acos(-1);
double mpi = 139.570/1000;//GeV
double mk = 493.677/1000;
double mp = 938.272/1000;
double mL = 1115.68/1000;
double mXi = 1321.71/1000;
double mXiStar = 1535/1000;
//	TFile* file = new TFile("run05000_DstTPCHelixTracking.root");
//	TTree* tree = (TTree*)file->Get("tpc");
/*
TFile* file = new TFile("SelectedHelixOld.root");
TTree* tree = (TTree*)file->Get("tree");
*/
void Clear(){
	ntTpc=0;
	isBeam->clear();
	pid->clear();
	charge->clear();
	chisqr->clear();
	mom0->clear();
	mom_vtx->clear();
	mom_vty->clear();
	mom_vtz->clear();
	closeDist->clear();
	vtx->clear();
	vty->clear();
	vtz->clear();
	helix_cx->clear();
	helix_cy->clear();
	helix_z0->clear();
	helix_r->clear();
	helix_dz->clear();
}
void SetHelixBranchAddress(TTree* tree);

void SetSummaryBranch(TTree* trout);
void Summary(){
	//	delete file;
	//	delete tree;
//	file = new TFile("run05000_DstTPCHelixTrackingOld.root");
	TFile* file = new TFile("run05000_DstTPCHelixTracking12.root");
	TTree* tree= (TTree*)file->Get("tpc");
	SetHelixBranchAddress(tree);

	int ent = tree->GetEntries();
	cout<<ent<<endl;
//	TFile* out = new TFile("SelectedHelixOld.root","recreate");
	TFile* out = new TFile("SelectedHelix12.root","recreate");
	TTree* trout = new TTree("tree","tree");
	SetSummaryBranch(trout);
	helixid=0;
	for(int i=0;i<ent;++i){
		Clear();
		tree->GetEntry(i);
		if(i%10000==0) cout<<i<<endl;
		if(ntTpc<1) continue;
		trout->Fill();
		helixid++;
	}
	out->Write();
}
void SetHelixBranchAddress(TTree* tree){
	
	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("evnum",&evnum);
	tree->SetBranchAddress("process_time",&process_time);
	tree->SetBranchAddress("ntTpc",&ntTpc);
	tree->SetBranchAddress("isBeam",&isBeam);
	tree->SetBranchAddress("chisqr",&chisqr);
	tree->SetBranchAddress("nclTpc",&nclTpc);
	tree->SetBranchAddress("cluster_de",&cluster_de);
	tree->SetBranchAddress("cluster_x",&cluster_x);
	tree->SetBranchAddress("cluster_y",&cluster_y);
	tree->SetBranchAddress("cluster_z",&cluster_z);
	tree->SetBranchAddress("cluster_size",&cluster_size);
	tree->SetBranchAddress("hough_flag",&hough_flag);
	tree->SetBranchAddress("helix_cx",&helix_cx);
	tree->SetBranchAddress("helix_cy",&helix_cy);
	tree->SetBranchAddress("helix_z0",&helix_z0);
	tree->SetBranchAddress("helix_r",&helix_r);
	tree->SetBranchAddress("helix_dz",&helix_dz);
	tree->SetBranchAddress("mom_vtx",&mom_vtx);
	tree->SetBranchAddress("mom_vty",&mom_vty);
	tree->SetBranchAddress("mom_vtz",&mom_vtz);
	tree->SetBranchAddress("vtx",&vtx);
	tree->SetBranchAddress("vty",&vty);
	tree->SetBranchAddress("vtz",&vtz);
	tree->SetBranchAddress("closeDist",&closeDist);
	tree->SetBranchAddress("mom0",&mom0);
	tree->SetBranchAddress("pid",&pid);
	tree->SetBranchAddress("charge",&charge);
	tree->SetBranchAddress("cluster_xf",&tcluster_x);
	tree->SetBranchAddress("cluster_yf",&tcluster_y);
	tree->SetBranchAddress("cluster_zf",&tcluster_z);
	tree->SetBranchAddress("track_cluster_x_center",&track_cluster_x_center);
	tree->SetBranchAddress("track_cluster_y_center",&track_cluster_y_center);
	tree->SetBranchAddress("track_cluster_z_center",&track_cluster_z_center);

}
void SetSummaryBranch(TTree* trout){
	trout->Branch("runnum",&runnum);
	trout->Branch("evnum",&evnum);
	trout->Branch("process_time",&process_time);
	trout->Branch("nclTpc",&nclTpc);
	trout->Branch("cluster_de",&cluster_de);
	trout->Branch("cluster_x",&cluster_x);
	trout->Branch("cluster_y",&cluster_y);
	trout->Branch("cluster_z",&cluster_z);
	trout->Branch("cluster_size",&cluster_size);
	trout->Branch("hough_flag",&hough_flag);
	trout->Branch("ntTpc",&ntTpc);
	trout->Branch("isBeam",&isBeam);
	trout->Branch("helix_cx",&helix_cx);
	trout->Branch("helix_cy",&helix_cy);
	trout->Branch("helix_z0",&helix_z0);
	trout->Branch("helix_r",&helix_r);
	trout->Branch("helix_dz",&helix_dz);
	trout->Branch("mom_vtx",&mom_vtx);
	trout->Branch("mom_vty",&mom_vty);
	trout->Branch("mom_vtz",&mom_vtz);
	trout->Branch("vtx",&vtx);
	trout->Branch("vty",&vty);
	trout->Branch("vtz",&vtz);
	trout->Branch("closeDist",&closeDist);
	trout->Branch("mom0",&mom0);
	trout->Branch("chisqr",&chisqr);
	trout->Branch("pid",&pid);
	trout->Branch("charge",&charge);
	trout->Branch("helixid",&helixid);
	trout->Branch("cluster_xf",&tcluster_x);
	trout->Branch("cluster_yf",&tcluster_y);
	trout->Branch("cluster_zf",&tcluster_z);
	trout->Branch("cluster_sizef",&tcluster_size);
	trout->Branch("track_cluster_x_center",&track_cluster_x_center);
	trout->Branch("track_cluster_y_center",&track_cluster_y_center);
	trout->Branch("track_cluster_z_center",&track_cluster_z_center);
}
