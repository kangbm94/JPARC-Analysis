#include "KKManager.hh"
void Acceptance(){
	TChain* chain = KM.GetPublicChain();
	TCanvas* c1 = new TCanvas("C1","C1",1200,600);
	c1->Divide(2,2);
	TH2D* Hist = new TH2D("PvsTheta","PvsTheta",100,0,30,100,0,3);
	TH1D* HistTh = new TH1D("Theta_PCut1","Theta_P<500",100,10,30);
	TH1D* HistTh2 = new TH1D("Theta_PCut2","Theta_P<530",100,10,30);
	TH1D* HistTh3 = new TH1D("Theta_PCut3","Theta_P<560",100,10,30);
	c1->cd(1);
	//	chain->Draw("pKurama[0]:thetaKurama[0]>>PvsTheta","qKurama[0]>0&&nKK==1&&inside[0]==1&&chisqrKurama[0]<20","col");
	chain->Draw("pKurama[0]:thetaKurama[0]>>PvsTheta","qKurama[0]>0&&nKK==1&&chisqrKurama[0]<200&&chisqrK18[0]<2&&chisqrK18[0]<200","col");
	c1->cd(2);
	chain->Draw("thetaKurama[0]>>Theta_PCut1","pKurama<0.5&&qKurama[0]>0&&nKK==1&&inside[0]==1&&chisqrKurama[0]<200","col");
	int peak = HistTh->GetMaximum();
	double ThetaCut = 14;
	TLine* ThetaCutV = new TLine(ThetaCut,0,ThetaCut,peak);
	ThetaCutV->SetLineWidth(2);
	ThetaCutV->SetLineColor(kRed);
	ThetaCutV->Draw("same");
	c1->cd(3);
	chain->Draw("thetaKurama[0]>>Theta_PCut2","pKurama<0.53&&qKurama[0]>0&&nKK==1&&inside[0]==1&&chisqrKurama[0]<200","col");
	peak = HistTh2->GetMaximum();
	ThetaCut = 13;
	TLine* ThetaCutV2 = new TLine(ThetaCut,0,ThetaCut,peak);
	ThetaCutV2->SetLineWidth(2);
	ThetaCutV2->SetLineColor(kRed);
	ThetaCutV2->Draw("same");
	c1->cd(4);
	chain->Draw("thetaKurama[0]>>Theta_PCut3","pKurama<0.56&&qKurama[0]>0&&nKK==1&&inside[0]==1&&chisqrKurama[0]<200","col");
	peak = HistTh2->GetMaximum();
	ThetaCut = 11;
	TLine* ThetaCutV3 = new TLine(ThetaCut,0,ThetaCut,peak);
	ThetaCutV3->SetLineWidth(2);
	ThetaCutV3->SetLineColor(kRed);
	ThetaCutV3->Draw("same");
}
void ViewScatterPlot(){
	TChain* chain = KM.GetPublicChain();
	TCanvas* c1 = new TCanvas("C1","C1",1200,600);
	c1->Divide(2,1);
	c1->cd(1);
	vector<double>target={4.2,0,19};
	vector<double>targetsize={20,10,30};
	TCut ChiCut="chisqrKurama[0]<50&&chisqrKurama[1]<50";
	//	TCut cut1= "nKK==2&&ntKurama==2&&abs(m2[0]-0.25)<0.1&&qKurama[0]<0&&pKurama[0]<1.4";
	TCut cut1= "nKK==2&&ntKurama==2&&abs(m2[1]-0.275)<0.125&&qKurama[1]<0&&pKurama[1]<1.4";
	//	TCut cut1= "nKK==2&&pKurama[1]<1.4";
	//	TCut cut2= "nKK==2&&ntKurama==2&&abs(m2[0]-0.25)<0.1&&qKurama[0]>0&&pKurama[0]<1.4";
	TCut cut2= "nKK==2&&ntKurama==2&&abs(m2[1]-0.275)<0.125&&qKurama[1]>0&&pKurama[1]<1.4";
	TCut vtc = VertexCut(target,targetsize,0)&&VertexCut(target,targetsize,1);
	TCut cd = "closeDist<4";
	//	vtc=vtc&&cd&&ChiCut;
	vtc="inside[0]==1&&inside[1]==1&&nKm==1";
	TH2D* h1 = new TH2D("KmTagged","KmTagged",100,-1,2,10,0,3);
	TH2D* h2 = new TH2D("KpTagged","KpTagged",100,-1,2,10,0,3);
	chain->Draw("pKurama[0]:sqrt(m2[0])*qKurama[0]>>KmTagged",cut1&&vtc,"colzgoff");
	//	chain->Draw("pKurama[1]:sqrt(m2[1])*qKurama[1]>>KmTagged",cut1&&vtc,"colzgoff");
	c1->cd(2);
	chain->Draw("pKurama[0]:sqrt(m2[0])*qKurama[0]>>KpTagged",cut2&&vtc,"colzgoff");
	//	chain->Draw("pKurama[1]:sqrt(m2[1])*qKurama[1]>>KpTagged",cut2&&vtc,"colzgoff");
	c1->cd(1);
	h1->Draw("col");
	c1->cd(2);
	h2->Draw("col");
}
void DrawMass(){
	TChain* chain = KM.GetPublicChain();
	int nb=1000;
	double m1= -1, m2 = 2;
	TH1D* hist = new TH1D("Mass","Mass",nb,m1,m2);
	TCut arg = "qKurama*sqrt(m2)>>Mass"; 
	TCut PCut = Form("pKurama<%f",1.4);
	TCut VertCut = Form("inside==1");
	TCut Cut = PCut&&VertCut;
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	chain->Draw(arg,Cut);
	double par[3];
	for(int i=0;i<6;i++){
		f_gauss[i] = new TF1(Form("f_gaus%d",i),"gaus");
	}
	for(int i=0;i<3;i++){
		GausMassFit(hist,i,-1,par);
		f_gauss[i]->SetParameters(par);
		cout<<par[1]<<endl;
		f_gauss[i]->Draw("same");
		f_gauss[i]->SetLineColor(kGreen);
		f_gauss[i]->SetRange(par[1]-0.05,par[1]+0.05);

		GausMassFit(hist,i,1,par);
		f_gauss[3+i]->SetParameters(par);
		f_gauss[3+i]->Draw("same");
		f_gauss[3+i]->SetRange(par[1]-0.05,par[1]+0.05);
	}

	/*
	double mean = Mass[0];
	double r1 = mean-MassWindow, r2 = mean+MassWindow;
	f_gaus->SetRange(r1,r2);
	f_gaus->SetParLimits(1,r1,r2);
	f_gaus->SetParLimits(2,0.01,0.1);
	hist->Fit("f_gaus","R");
	
	mean = -Mass[0];
	r1 = mean-MassWindow, r2 = mean+MassWindow;
	f_gaus->SetRange(r1,r2);
	f_gaus->SetParLimits(1,r1,r2);
	f_gaus->SetParLimits(2,0.01,0.1);
	
	hist->Fit("f_gaus","R");
*/

}
