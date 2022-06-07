double dt[10][500];
int tdc[10][500];
double tot[10][500];
int wire[10][500];
int nhit[10];
double theta[10];
bool TrackFlag[10][500];
TTree* gtree;
TChain* gchain;
TTree* ghist;
TFile* gfile;
TFile* gHistfile;
TTree* gpartree;
TString gFilename;
double T0RangeMin[10]={420,420,420,420,420,420,930,930,930,920};
double T0RangeMax[10]={467,467,466,465,469,468,975,975,975,975};
double dlpitch[10] = {3,3,3,3,3,3,5,5,5,5};


double chisqr[10];
int nw,nb;
double p0[70];double p1[70];double p2[70];double p0err[70];double p1err[70];double p2err[70];double chi2[70];int ndf[70];
double d_p[6];
void InitParameters(){
	nw=0;
	nb=0;
	for(int i=0;i<70;i++){
		p0[i]=0;
		p1[i]=0;
		p2[i]=0;
		p0err[i]=0;
		p1err[i]=0;
		p2err[i]=0;
		chi2[i]=0;
		ndf[i]=0;
	}
}


double FitFunc(double* x, double* p){
//	return (p[0]*(1+erf(p[1]*(x[0]-p[2]) )  )/2+p[3])*Step(p[4]-x[0]);
//	return (p[0]*(1+erf(p[1]*(x[0]-p[2]) )  )/2+p[3]);
	return (p[0]*(1+erf(p[1]*(x[0]-p[2]) )  )/2);
}
TF1* func = new TF1("func","FitFunc",0,1500,3);


static double lsb=0.833;
TTree* AssignTree(TString filename){
	TTree* tree;
	TFile* file = new TFile(filename,"READ");
	tree = (TTree*)file->Get("sdcin");
	tree->SetBranchAddress("dt",dt);
	tree->SetBranchAddress("tdc",tdc);
	tree->SetBranchAddress("tot",tot);
	tree->SetBranchAddress("wire",wire);
	tree->SetBranchAddress("TrackFlag",TrackFlag);
	tree->SetBranchAddress("nhit",nhit);
	tree->SetBranchAddress("theta",theta);
	tree->SetBranchAddress("chisqr",chisqr);
	return tree;
}
TChain* AssignChain(TString filename){
	TChain* chain=new TChain("sdcin");
//	TFile* file = new TFile(filename,"READ");
	chain->Add(filename);
	chain->SetBranchAddress("dt",dt);
	chain->SetBranchAddress("tdc",tdc);
	chain->SetBranchAddress("tot",tot);
	chain->SetBranchAddress("wire",wire);
	chain->SetBranchAddress("TrackFlag",TrackFlag);
	chain->SetBranchAddress("nhit",nhit);
	chain->SetBranchAddress("theta",theta);
	chain->SetBranchAddress("chisqr",chisqr);
	return chain;
}
int Time_to_TDC(double time,double p0){
 return int(p0-(time)/(lsb));
}
double TDC_to_Time(int tdc,double p0){
 return lsb*(p0-tdc);
}
double SDCRise(double* x,double* p){
	double val;
	return p[0]+p[1]*tanh( (x[0]-p[2])/p[3]);
}
double pol5s(double* x,double* p){
	double val=0;
	for(int i=0;i<5;i++){
		val+=p[i+1]*pow(x[0]-p[0],i+1);
	}
	return val;
}
double pol5_chev(double* x,double* p){
	double val=0;
	for(int i=0;i<6;i++){
		if(i<4){
			val+=p[6]*p[i]*chebyshev(x[0]/p[5],i);
		}
		else if(i==4){
			val+=p[6]*(p[2]-p[0])*chebyshev(x[0]/p[5],i);
		}
		else {
			val+=p[6]*p[i-1]*chebyshev(x[0]/p[5],i);
		}
	}
	return val;
}
TF1* pol5f = new TF1("pol5f","pol5s",-50,250,6);
TF1* pol5_chevf = new TF1("pol5_chevf","pol5_chev",0,120,7);
TF1* SDCRisef= new TF1("SDCRisef","SDCRise",0,120,4);
void LoadTree(TString filename){
	gtree= AssignTree(filename);
	gFilename=filename;
}
void LoadChain(TString filename){
	gchain= AssignChain(filename);
	gFilename=filename;
}
void LoadFile(TString filename){
	gfile=new TFile("Hist"+filename,"RECREATE");
}
void MakeMTHist(int layer){
	gfile->cd();
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	int l1 = 2*((layer-1)/2),l2=l1+1;
	double sdc1t1=250,sdc1t2=500;
	double sdc2t1=700,sdc2t2=1100;
	int sdc1nw=	64;
	int sdc2nwx=	70;
	int sdc2nwy=	40;
	int nbin1=(sdc1t2-sdc1t1);
	int nbin2=(sdc2t2-sdc2t1);
	double t1=300,t2=0;
	int nw,nb;
	double theta_cut=12.,chicut=1.5;
	double totcut1[10]={30,30,30,30,30,30,130,150,180,130};//noise cut
	double totcut2[10]={90,90,90,90,90,90,250,270,270,250};//noise cut
	if(layer<7){	t2=sdc1t2;nb=nbin1;t1=sdc1t1;	nw=sdc1nw;}
	else if(layer<9){	t2=sdc2t2;nb=nbin2;t1=sdc2t1;	nw=sdc2nwx;}
	else{	t2=sdc2t2;nb=nbin2;t1=sdc2t1;	nw=sdc2nwy;}
	vector<TH1D*> tdchist(nw);
	vector<TString> ht(nw);
	TString htt= Form("MeanTimeL%dL%d",l1,l2);
	TH1D* HistMT= new TH1D(htt,htt,nb*2,t1,t2);
	TString drw = Form("tdc[%d]/2+tdc[%d]/2>>",l1,l2)+htt;
	cout<<drw<<endl;
	TCut cuts = "ntrack==1";
	gchain->Draw(drw,cuts);
	double r1=T0RangeMin[l1]-40,r2=T0RangeMin[l1]+10;
	/*
	func_dgaus->SetRange(r1,r2);
	func_dgaus->SetParLimits(0,r1-10,T0RangeMax[l1]);
	func_dgaus->SetParLimits(1,2,20);
	func_dgaus->SetParLimits(3,-30,-10);
	func_dgaus->SetParLimits(4,2,20);
	*/
//	HistMT->Fit("func_dgaus","R");
}

void MakeMTHist(){
	gfile=new TFile("HistMT"+gFilename,"RECREATE");
	for(int i=1;i<6;i++){
		MakeMTHist(i*2);
	}
	gfile->Write();
}



void SdcIn(){
	//05047,05330,05340,05341,05342,05343
	/*
	LoadChain("SdcInParameter05340.root");
	gchain->Add("SdcInParameter05341.root");
	gchain->Add("SdcInParameter05342.root");
	gchain->Add("SdcInParameter05343.root");
	*/
	LoadChain("SdcInParameter_TrackSelection_05047_0.root");
	//	LoadTree("SdcInParameter05395.root");
	cout<<"MakeT0Hist()"<<endl;
	cout<<"MakeT0Hist()"<<endl;	
	cout<<"FitT0()"<<endl;
	cout<<"IntegrateDriftTime()"<<endl;
}
void IntegrateDriftTime(int layer){
	int l1 = 2*((layer-1)/2),l2=l1+1;
	int cnt=0;
	double sdc1t1=-30,sdc1t2=180;
	double sdc2t1=-50,sdc2t2=250;
	int sdc1nw=	64;
	int sdc2nwx=	70;
	int sdc2nwy=	40;
	int nbin1=(sdc1t2-sdc1t1)*6/5;//6/5 = 1/0.833 = 1/lsb
	int nbin2=(sdc2t2-sdc2t1)*6/5;
	double t1=300,t2=0;
	int nw,nb;
	double theta_cut=12.,chicut=1.5;
	double mtcut[5]={409.96,408.98,411.7,890,800};
//	TCut cut = Form("(tdc[%d]+tdc[%d])/2>%f",l1,l2,mtcut[int((layer-1)/2)]);
	TCut cut ="nlayer==10";
	//Form("(tdc[%d]+tdc[%d])/2>%f",l1,l2,mtcut[int((layer-1)/2)]);
	cout<<cut<<endl;
	if(layer<7){	t2=sdc1t2;nb=nbin1;t1=sdc1t1;	nw=sdc1nw;}
	else if(layer<9){	t2=sdc2t2;nb=nbin2;t1=sdc2t1;	nw=sdc2nwx;}
	else{	t2=sdc2t2;nb=nbin2;t1=sdc2t1;	nw=sdc2nwy;}
	TString dtht = Form("DT_L%d",layer);
	TH1D* dthist = new TH1D(dtht,dtht,nb,t1,t2);
	TString drw = Form("dt[%d]>>",layer-1)+dtht;
	cout<<drw<<endl;
	gchain->Draw(drw,cut);	
}

void IntegrateDriftTime(){
	gfile = new TFile("SdcInDTHist_05047.root","recreate");
	TCanvas* c1 = new TCanvas("c1","c1",1500,800);
	c1->Divide(4,3);
	for(int i=1;i<=10;i++){
		cout<<Form("Layer %d",i)<<endl;
		c1->cd(i);
		IntegrateDriftTime(i);
	}
	gfile->Write();
}
void MakeDLDTHist(int layer){
//	TCanvas* c1 = new TCanvas(Form("DLDT%d",layer),Form("DLDT%d",layer),1200,600);
	//	gfile = new TFile("SdcInDTHist_05340.root","READ");
	TH1D* h=(TH1D*)gfile->Get(Form("DT_L%d",layer));
	int PP = h->GetMaximumBin();
	int Peak = h->GetBinContent(PP);
	int nb = h->GetNbinsX();
	int HP=Peak/4;//2 -> Wider
	int End=0;
	int Start=0;
	bool SFlg=true;
	int buf[1000];	
	double bw,bs,be;//binwidth,startbin, endbin
	for(int i=nb;i>0;i--){
		if(h->GetBinContent(i)>HP){
			End=i;
			bw=h->GetBinWidth(i);
			be=h->GetBinCenter(i)+bw/2;
			i=0;
		}
	}
	buf[0]=h->GetBinContent(1);
	buf[1]=h->GetBinContent(2);
	buf[2]=h->GetBinContent(3);
	for(int i=3;i<nb;i++){
		SFlg=true;
		buf[i]=h->GetBinContent(i);
		for(int j=0;j<3;j++){
//			if(buf[i]<=buf[i-1-j]){
			if(buf[i]<1){
				SFlg=SFlg*false;
			}
		}
		if(SFlg){
			Start=i;
			bs=h->GetBinCenter(i)-bw/2;
			cout<<"Starting Time: "<<h->GetBinCenter(Start)<<endl;
			cout<<"Ending Time: "<<h->GetBinCenter(End)<<endl;
			i=nb;
			continue;
		}
	}
	int nbin = End-Start;
	cout<<Form("Layer %d nbins: %d",layer,nbin)<<endl;
	int dli[1000];
	double dls[1000];double dle[1000];double dts[1000];double dte[1000];
	TString ht = Form("DLDT_L%d",layer);
	pol5f->SetRange(bs,be);
	pol5f->SetParLimits(0,-20,0);
	pol5f->SetParLimits(1,-1e-1,1e-1);
	pol5f->SetParLimits(2,-1e-2,1e-2);
	pol5f->SetParLimits(3,-1e-4,1e-4);
	pol5f->SetParLimits(4,-1e-6,1e-6);
	pol5f->SetParLimits(5,-1e-8,1e-8);
//	TH1D* hist_dldt = new TH1D(ht,ht,nbin,bs,be);
	int ents=0;
	for(int i=0;i<nbin;i++){
			ents+=h->GetBinContent(Start+i);
	}
	for(int i=0;i<nbin;i++){
		if(i==0){
			dli[i]=h->GetBinContent(Start+i);
		}else{
			dli[i]=dli[i-1]+h->GetBinContent(Start+i);}
			dls[i]=dlpitch[layer-1]*(dli[i])/ents;
			dle[i]=dlpitch[layer-1]/nbin;
			dts[i]=h->GetBinCenter(Start+i);
			dte[i]=h->GetBinWidth(Start+i)/sqrt(12);
	}
	cout<<"Normalization: "<<dli[nbin-1]/ents<<endl;
	TGraphErrors* ge = new TGraphErrors(nbin,dts,dls,dte,dle);
	gHistfile->cd();
	ge->SetTitle(Form("DLDT_%d",layer));
	pol5f->SetLineWidth(6);
	ge->Draw("AP");
	ge->Fit("pol5f","QR");
	ge->Fit("pol5f","R");
	for(int i=0;i<6;i++){
		d_p[i]=pol5f->GetParameter(i);
	}
}

void MakeDLDTHist(){
	gfile = new TFile("SdcInDTHist_05047.root","READ");
	gHistfile = new TFile("SdcInDLDT_05047.root","RECREATE");	
	TCanvas* c1 = new TCanvas("c1","c1",1500,800);
	c1->Divide(4,3);
	fstream f;
	string ftitle = "DCDRFT_05047";
	f.open(ftitle,fstream::out);
	for(int i=1;i<=10;i++){
		c1->cd(i);
		cout<<Form("%dth layer",i)<<endl;
		MakeDLDTHist(i);
		if(i==1){
			f<<"###	SDC1"<<endl;
		}
		if(i==7){
			f<<"###	SDC2"<<endl;
		}
		f<<i<<"	"<<0<<"	"<<6<<" "<<6;
		for(int j=0;j<6;j++){
			f<<" "<<d_p[j];
		}
		f<<endl;
	}
	gHistfile->Write();
}
void FitT0(int layer){
	InitParameters();
	int cnt=0;
	double sdc1t1=250,sdc1t2=500;
	double sdc2t1=700,sdc2t2=1100;
	int sdc1nw=	64;
	int sdc2nwx=	70;
	int sdc2nwy=	40;
	int nbin1=(sdc1t2-sdc1t1);
	int nbin2=(sdc2t2-sdc2t1);
	double t1=300,t2=0;
	double p2ErrCut=1.99*lsb;
	if(layer<7){	t2=sdc1t2;nb=nbin1;t1=sdc1t1;	nw=sdc1nw;}
	else if(layer<9){	t2=sdc2t2;nb=nbin2;t1=sdc2t1;	nw=sdc2nwx;}
	else{	t2=sdc2t2;nb=nbin2;t1=sdc2t1;	nw=sdc2nwy;}
	vector<TH1D*> tdchist(nw);
	vector<TString> ht(nw);
	double wire[70];double we[70];
	for(int i=0;i<nw;i++){
		ht[i]=Form("L%dW%d",layer,i+1);
		tdchist[i]= (TH1D*)gHistfile->Get(ht[i]);
		wire[i]=i+1;
		we[i]=0;
	}
//	cout<<tdchist[i]->GetTitle()<<endl;
	int niter=1;
	double p1Fixed=0,p2Mean=0;
	int p1Cnt=0;
	double rm = T0RangeMin[layer-1],rM=T0RangeMax[layer-1];
	TCanvas* c1= new TCanvas(Form("L%d",layer),Form("L%d",layer),1800,1200);
	if(layer<9){
		c1->Divide(10,7);
	}
	else{
		c1->Divide(8,5);
	}
	for(int i=0;i<nw;i++){
		c1->cd(i+1);
		for(int j=0;j<niter;j++){
			func->SetRange(rm,rM);
			int PP = tdchist[i]->FindBin(rm);//PeakPostion
			int TP = tdchist[i]->FindBin(rM);//TailPostion
			int peak = tdchist[i]->GetBinContent(PP);
			int tail = tdchist[i]->GetBinContent(TP);
			cout<<Form("Layer: %d, wire: %d,  Peak: %d, Tail: %d",layer,i+1,peak,tail)<<endl;
			cout<<Form("Range: (%d,%d)",int(rm),int(rM))<<endl;
			func->SetParLimits(0,0.8*peak,1.5*peak);
			func->SetParLimits(1,-0.12,-0.09);
			func->SetParLimits(2,rm,rM);
			func->SetParLimits(3,0,0.3*peak);
//			func->FixParameter(3,0);
//			func->SetParLimits(4,rM-5,rM+5);
			tdchist[i]->Fit("func","QR");
			p0[i]=func->GetParameter(0);
			p0err[i]=func->GetParError(0);
			p1[i]=func->GetParameter(1);
			p1err[i]=func->GetParError(1);
			p2[i]=func->GetParameter(2);
			p2err[i]=func->GetParError(2);
			ndf[i]=func->GetNDF();
			chi2[i]=p2[i]-1/p1[i];
			cout<<p2[i]<<endl;
			if(p2err[i]<p2ErrCut){
				p1Fixed+=p1[i];
				p2Mean+=p2[i];
				p1Cnt++;
			}
		}
	}
	p1Fixed=p1Fixed/p1Cnt;p2Mean=p2Mean/p1Cnt;
//	TCanvas* c2= new TCanvas(Form("L2%d",layer),Form("L2%d",layer),1800,1200);
	for(int i=0;i<nw;i++){
		c1->cd(i+1);
		for(int j=0;j<niter;j++){
			int PP = tdchist[i]->FindBin(rm);//PeakPostion
			int TP = tdchist[i]->FindBin(rM);//TailPostion
			int peak = tdchist[i]->GetBinContent(PP);
			int tail = tdchist[i]->GetBinContent(TP);
//			cout<<Form("Wire%d, peak: %f",i,PL)<<endl;
//			func->SetRange(p2Mean-10,p2Mean+8);
			func->SetParLimits(2,rm,rM);
			func->SetParLimits(0,0.8*peak,1.5*peak);
			func->FixParameter(1,p1Fixed);
//			func->SetParLimits(1,p1Fixed*1.2,p1Fixed*0.8);
			func->SetParLimits(2,p2Mean-5,p2Mean+5);
			func->SetParLimits(3,0,0.3*peak);
//			func->FixParameter(3,0);
			func->SetParLimits(4,rM-5,rM+5);
			tdchist[i]->Fit("func","QR");
//			tdchist[i]->Fit("gaus","Q");
		
			p0[i]=func->GetParameter(0);
			p0err[i]=func->GetParError(0);
			p1[i]=func->GetParameter(1);
			p1err[i]=func->GetParError(1);
			p2[i]=func->GetParameter(2);
			p2err[i]=func->GetParError(2);
			ndf[i]=func->GetNDF();
			chi2[i]=p2[i]-1/p1[i];
//			chi2[i]=(func->GetChisquare())/ndf[i];
			}
	}
	TCanvas* c2= new TCanvas(Form("Par%d",layer),Form("Par%d",layer),1800,1200);
	c2->Divide(2,2);
	TGraphErrors* g[4];
	g[0]=new TGraphErrors(nw,wire,p0,we,p0err);
	g[1]=new TGraphErrors(nw,wire,p1,we,p1err);
	g[2]=new TGraphErrors(nw,wire,p2,we,p2err);
	g[3]=new TGraphErrors(nw,wire,chi2,we,we);
	g[0]->SetTitle("p0");
	g[1]->SetTitle("p1");
	g[2]->SetTitle("p2");
	g[3]->SetTitle("chi2");
	for(int i=0;i<4;i++){
		c2->cd(i+1);
		g[i]->Draw();
	}
}

void FitT0(){
	fstream f;
	string ftitle = "DCTdcCalib_05047";//rn
	f.open(ftitle,fstream::out);
	TFile* pars = new TFile("T0"+gFilename,"RECREATE");
	gpartree=new TTree("Parameters","tree");
	gpartree->Branch("nw",&nw,"nw/I");
	gpartree->Branch("p0",p0,"p0[nw]/D");
	gpartree->Branch("p1",p1,"p1[nw]/D");
	gpartree->Branch("p2",p2,"p2[nw]/D");
	gpartree->Branch("p0err",p0err,"p0err[nw]/D");
	gpartree->Branch("p1err",p1err,"p1err[nw]/D");
	gpartree->Branch("p2err",p2err,"p2err[nw]/D");
	gpartree->Branch("chi2",chi2,"chi2[nw]/D");
	gpartree->Branch("ndf",ndf,"ndf[nw]/I");
	gHistfile = new TFile("Hist"+gFilename,"READ");
	string lt[10]={"### SDC1-V1","### SDC1-V2","### SDC1-X1","### SDC1-X2","### SDC1-U1","### SDC1-U2","### SDC2-X1","### SDC2-X2","### SDC2-Y1","### SDC2-Y2"};

	for(int i=1;i<=10;i++){
		f<<lt[i-1]<<endl;
		f<<"#"<<endl;
		FitT0(i);
		for(int j=0;j<nw;j++){
			f<<i<<"	"<<j+1<<"	"<<0.833<<" "<<-p2[j]<<endl;
		}
		gpartree->Fill();
	}
	pars->Write();
}
void MakeT0Hist(int layer){
	cout<<Form("Processing Layer %d",layer)<<endl;
	int l1 = 2*((layer-1)/2),l2=l1+1;
	int cnt=0;
	double sdc1t1=250,sdc1t2=500;
	double sdc2t1=700,sdc2t2=1100;
	int sdc1nw=	64;
	int sdc2nwx=	70;
	int sdc2nwy=	40;
	int nbin1=(sdc1t2-sdc1t1);
	int nbin2=(sdc2t2-sdc2t1);
	double t1=300,t2=0;
	int nw,nb;
	double theta_cut=12.,chicut=1.5;
	double totcut1[10]={30,30,30,30,30,30,130,150,180,130};//noise cut
	double totcut2[10]={90,90,90,90,90,90,250,270,270,250};//noise cut
//	double mtcut[5]={409.96,408.98,411.7,890,800};
	if(layer<7){	t2=sdc1t2;nb=nbin1;t1=sdc1t1;	nw=sdc1nw;}
	else if(layer<9){	t2=sdc2t2;nb=nbin2;t1=sdc2t1;	nw=sdc2nwx;}
	else{	t2=sdc2t2;nb=nbin2;t1=sdc2t1;	nw=sdc2nwy;}
	vector<TH1D*> tdchist(nw);
	vector<TString> ht(nw);
	
	for(int i=0;i<nw;i++){
		ht[i]=Form("L%dW%d",layer,i+1);
		tdchist[i]= new TH1D(ht[i],ht[i],nb+1,t1,t2);
	}
	TObjString s(Form("Par%d",layer));
  gfile->WriteObject(&s,Form("Layer%d_Cuts=(theta<%f,chisqr<%f,%d<tot<%d)",layer,theta_cut,chicut,int(totcut1[layer]),int(totcut2[layer])));
	int ent = gchain->GetEntries();
	cout<<Form("Layer %d",layer)<<endl;
//	ent=ent/1000;
	for(int i=0;i<ent;i++){
		gchain->GetEntry(i);
		Indicator(i,ent);
		for(int j=0;j<nhit[layer-1];j++){
//			if(TrackFlag[0][j]==0)continue;
//			if((tdc[l1][j]+tdc[l2][j])/2<mtcut[int((layer-1)/2)]){continue;}
			int w = wire[layer-1][j]-1;
			if(w<0||w>nw)continue;
			tdchist[w]->Fill(tdc[layer-1][j]);
			cnt++;
		}
	}
	cout<<Form("Layer%d: Accepted %d",layer,cnt)<<endl;
}

void MakeT0Hist(){
	LoadFile(gFilename);
	for(int i=1;i<=10;i++){
		MakeT0Hist(i);
	}
	gfile->Write();
}

