
TVector3
VertexPoint(const TVector3& Xin, const TVector3& Xout,
            const TVector3& Pin, const TVector3& Pout)
{
  Double_t xi=Xin.x(), yi=Xin.y(), xo=Xout.x(), yo=Xout.y();
  Double_t ui=Pin.x()/Pin.z(), vi=Pin.y()/Pin.z();
  Double_t uo=Pout.x()/Pout.z(), vo=Pout.y()/Pout.z();

  Double_t z=((xi-xo)*(uo-ui)+(yi-yo)*(vo-vi))/
    ((uo-ui)*(uo-ui)+(vo-vi)*(vo-vi));
  Double_t x1=xi+ui*z, y1=yi+vi*z;
  Double_t x2=xo+uo*z, y2=yo+vo*z;
  Double_t x = 0.5*(x1+x2);
  Double_t y = 0.5*(y1+y2);
//  if(std::isnan(x) || std::isnan(y) || std::isnan(z))
 //   return TVector3(nan, nan, nan);

  return TVector3(x, y, z);

}

	Double_t
CloseDist(const TVector3& Xin, const TVector3& Xout,
          const TVector3& Pin, const TVector3& Pout)
{
  Double_t xi=Xin.x(), yi=Xin.y(), xo=Xout.x(), yo=Xout.y();
  Double_t ui=Pin.x()/Pin.z(), vi=Pin.y()/Pin.z();
  Double_t uo=Pout.x()/Pout.z(), vo=Pout.y()/Pout.z();

  Double_t z=((xi-xo)*(uo-ui)+(yi-yo)*(vo-vi))/
    ((uo-ui)*(uo-ui)+(vo-vi)*(vo-vi));
  Double_t x1=xi+ui*z, y1=yi+vi*z;
  Double_t x2=xo+uo*z, y2=yo+vo*z;

  return TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}


void CheckUndetected(){
	TFile* file = new TFile("undetected.root");
	TTree* tree = (TTree*)file->Get("tree");
	int evn;
	int runnum = 5641;
	tree->SetBranchAddress("evnum",&evn);
	TFile* filek = new TFile(Form("run0%d_DstHSKKAna.root",runnum));
	TTree* treek = (TTree*)filek->Get("kk");
	int nKm,nKp,nKK;
	double KMPX[5],KMPY[5],KMPZ[5];
	double KPPX[5],KPPY[5],KPPZ[5];
	double xkm[5],ykm[5],xkp[5],ykp[5];
	double chisqrKurama[5],closeDist[5],MissMass[5],m2[5],pKurama[5],qKurama[5];
	double chisqrK18[5];
	int inside[5];
	treek->SetBranchAddress("nKm",&nKm);
	treek->SetBranchAddress("nKp",&nKp);
	treek->SetBranchAddress("nKK",&nKK);
	treek->SetBranchAddress("m2",m2);
	treek->SetBranchAddress("inside",inside);
	treek->SetBranchAddress("chisqrKurama",chisqrKurama);
	treek->SetBranchAddress("chisqrK18",chisqrK18);
	treek->SetBranchAddress("pKurama",pKurama);
	treek->SetBranchAddress("qKurama",qKurama);
	treek->SetBranchAddress("closeDist",closeDist);
	treek->SetBranchAddress("KMPX",KMPX);
	treek->SetBranchAddress("KMPY",KMPY);
	treek->SetBranchAddress("KMPZ",KMPZ);
	treek->SetBranchAddress("KPPX",KPPX);
	treek->SetBranchAddress("KPPY",KPPY);
	treek->SetBranchAddress("KPPZ",KPPZ);
	treek->SetBranchAddress("xkm",xkm);
	treek->SetBranchAddress("ykm",ykm);
	treek->SetBranchAddress("xkp",xkp);
	treek->SetBranchAddress("ykp",ykp);
	treek->SetBranchAddress("MissMassCorr",MissMass);
	TH2D* hvert = new TH2D("vert","vert",50,-200,400,50,-90,90);
	TH1D* hcd = new TH1D("cd","cd",100,0,100);
	TH1D* hchi = new TH1D("chi","chi",1000,00,1000);
	TH1D* hm = new TH1D("Missmass","Missmass",100,1,2);
	TH2D* hm2 = new TH2D("pKurama:mass","pKurama:qKurama*mass",100,-1,3,100,0,2);
	TH1I* hflag = new TH1I("Flag","Flag",10,0,10);

	int ent = tree->GetEntries();
	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
		treek->GetEntry(evn);
		TVector3 XKM,PKM,XKP,PKP;
		XKM = TVector3(xkm[0],ykm[0],0);
		PKM = TVector3(KMPX[0],KMPY[0],KMPZ[0]);
		XKP = TVector3(xkp[0],ykp[0],0);
		PKP = TVector3(KPPX[0],KPPY[0],KPPZ[0]);
		auto Vert = VertexPoint(XKM,XKP,PKM,PKP);
//		double cd = CloseDist(XKM,XKP,PKM,PKP);
		double x = Vert.x();
		double y = Vert.y();
		double z = Vert.z();
//		cout<<Form("Vertex (%f,%f,%f)",x,y,z)<<endl;
		if(chisqrK18[0]>20)cout<<"Warning"<<endl;
		if(nKK!=1){
			cout<<"Event : "<<evn<<endl;
			cout<<"nKK = "<<nKK<<endl;
//			continue;
			hflag->Fill(0);
		}
		if(inside[0]!=1){
			cout<<"Event : "<<evn<<endl;
			cout<<"inside = "<<inside[0]<<endl;
//			continue;
			hflag->Fill(1);
		}
		if(chisqrKurama[0]>200){
			cout<<"Event : "<<evn<<endl;
			cout<<"chisqrKurama = "<<chisqrKurama[0]<<endl;
//			continue;
			hflag->Fill(2);
		}
		if(pKurama[0]>1.4){
			cout<<"Event : "<<evn<<endl;
			cout<<"pKurama = "<<pKurama[0]<<endl;
//			continue;
			hflag->Fill(3);
		}
		hvert->Fill(z,x);
		hcd->Fill(closeDist[0]);
		hchi->Fill(chisqrKurama[0]);
		hm2->Fill(qKurama[0]*sqrt(m2[0]),pKurama[0]);
		hm->Fill(MissMass[0]);
	}
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(3,2);
	c1->cd(1);
	hvert->Draw("colz");
	cout<<hvert->GetEffectiveEntries()<<endl;
	c1->cd(2);
	hcd->Draw("colz");
	c1->cd(3);
	hchi->Draw();
	c1->cd(4);
	hm2->Draw("colz");
	c1->cd(5);
	hm->Draw();
	c1->cd(6);
	hflag->Draw();
}	
