void Compare(){
	gStyle->SetOptStat(0);
	TFile* kbmf = new TFile("run05641_XiEvents_KBM.root");
//	TFile* jwsf = new TFile("run05641_XiEvents_JWS.root");
	TFile* jwsf = new TFile("ws.root");
	TFile* cswf = new TFile("run05641_XiEvents_CSW.root");
	TTree* kbmt = (TTree*)kbmf->Get("tree");
	TTree* jwst = (TTree*)jwsf->Get("tree");
	TTree* cswt = (TTree*)cswf->Get("tree");
	int evk,evj,evc;
	kbmt->SetBranchAddress("evnum",&evk);
//	jwst->SetBranchAddress("evnum",&evj);
	jwst->SetBranchAddress("evtid",&evj);
//	jwst->SetBranchAddress("evnum",&evj);
	cswt->SetBranchAddress("evnum",&evc);
	vector<int>kj;
	vector<int>jc;
	vector<int>ck;
	vector<int>kjc;
	TH1I* hist = new TH1I("evs","evs",8,0,8);
	TFile* File = new TFile("undetected.root","recreate");
	TTree* tree = new TTree("tree","tree");
	tree->Branch("evnum",&evk);
	for(int i=0;i<kbmt->GetEntries();++i){
		kbmt->GetEntry(i);
		for(int j=0;j<jwst->GetEntries();++j){
			jwst->GetEntry(j);
			if(evk == evj){
				kj.push_back(evk);	
				break;
			}
		}	
	}
	for(int i=0;i<jwst->GetEntries();++i){
		jwst->GetEntry(i);
		for(int j=0;j<cswt->GetEntries();++j){
			cswt->GetEntry(j);
			if(evj == evc){
				jc.push_back(evj);	
				break;
			}
		}	
	}
	for(int i=0;i<kbmt->GetEntries();++i){
		kbmt->GetEntry(i);
		for(int j=0;j<cswt->GetEntries();++j){
			cswt->GetEntry(j);
			if(evk == evc){
				ck.push_back(evk);	
				break;
			}
		}	
	}
	for(auto ev:kj){
		for(int j=0;j<cswt->GetEntries();++j){
			cswt->GetEntry(j);
			if(ev == evc){
				kjc.push_back(ev);	
				break;
			}
		}	
	}
	int nk = kbmt->GetEntries();
	int nj = jwst->GetEntries();
	int nc = cswt->GetEntries();
	int nkj = kj.size();
	int njc = jc.size();
	int nck = ck.size();
	int nkjc = kjc.size();
	
	for(int i=0;i<nj;++i){
		jwst->GetEntry(i);
		bool acpt = true;
		for(int j = 0;j<nk;++j){
			kbmt->GetEntry(j);
			if(evj == evk){
				acpt = false;
				break;
			}
		}
		evk = evj;
//		if(acpt)tree->Fill();	
	}
	for(int i=0;i<nc;++i){
		cswt->GetEntry(i);
		bool acpt = true;
		for(auto ev : jc){
//			if(evc == ev)acpt = false;
		}
		for(int j = 0;j<nk;++j){
			kbmt->GetEntry(j);
			if(evc == evk){
				acpt = false;
				break;
			}
		}
		evk = evc;
		if(acpt)tree->Fill();	
	}
	File->Write();
	hist->SetBinContent(1,nk-nkj-nck+nkjc);
	hist->SetBinContent(2,nj-nkj-njc+nkjc);
	hist->SetBinContent(3,nc-njc-nck+nkjc);
	hist->SetBinContent(4,nkj-nkjc);
	hist->SetBinContent(5,njc-nkjc);
	hist->SetBinContent(6,nck-nkjc);
	hist->SetBinContent(7,nkjc);
	hist->SetBinContent(8,nk+nj+nc-nkj-nck-njc+nkjc);
	hist->Draw();
	cout<<Form("bm,ws,sw = (%d,%d,%d)",nk,nj,nc)<<endl;
	TLegend* L = new TLegend(0.1,0.6,0.6,0.9);
	L->Draw();
	L->AddEntry((TObject*)0,"0,1,2 = KBM/JWS/CSW Only");
	L->AddEntry((TObject*)0,"3,4,5 = K&J/J&C/C&K Only");
	L->AddEntry((TObject*)0,"6 = All AND");
	L->AddEntry((TObject*)0,"7 = All OR");
	for(int i=0;i<8;++i){
		cout<<Form("ent%d = %f",i,hist->GetBinContent(i+1))<<endl;
	}
}
