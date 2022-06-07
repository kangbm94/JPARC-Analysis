#include "EasirocManager.hh"
EasirocManager EM;
void Easiroc(){
	EM.LoadFile("./rootfiles/CH2/BigEasiroc.root");
}
void BftT0(){
	double T_min=650,T_max=700;
	
	TH1* h[2*bftnseg];
	TLine* Line[2*bftnseg];	
	vector<int> ID = {110,0,0,1,0};//Cid, PlID, SegID, AorT, UorD. For BFT, Set PlID instead of UorD
	vector<double> Param = {0,-1};

	EM.MakeParameterFile("./Params/BFT_T0.txt");
	EM.WriteComment("#########################################################");
	EM.WriteComment("# BFT-U Tdc                                             #");
	EM.WriteComment("#########################################################");
	ID[1]=0;
	TCanvas* Canv_Up = new TCanvas("UP","UP",1600,1200);
	Canv_Up->Divide(16,10);

	for(int i=0;i<bftnseg;i++){
		ID[2]=i;
		h[i]= (TH1*)EM.DrawBftTDC(0,i+1);
		h[i]->SetAxisRange(	T_min,T_max);
		double PeakPosition = h[i]->GetBinCenter(h[i]->GetMaximumBin());
		double PeakHeight = h[i]->GetMaximum();
		Line[i]= new TLine(PeakPosition,0,PeakPosition,PeakHeight);
		Line[i]->SetLineWidth(2);
		Line[i]->SetLineColor(kBlue);
		Canv_Up->cd(i+1);
		Param[0]=PeakPosition;
		EM.WriteParameter(ID,Param);
		h[i]->Draw();
		Line[i]->Draw("SAME");
	}
	EM.WriteComment("#########################################################");
	EM.WriteComment("# BFT-D Tdc                                             #");
	EM.WriteComment("#########################################################");
	ID[1]=1;
	TCanvas* Canv_Dn = new TCanvas("Dn","Dn",1600,1200);
	
	Canv_Dn->Divide(16,10);
	for(int i=0;i<bftnseg;i++){
		ID[2]=i;
		cout<<i<<endl;
		h[bftnseg+i]= (TH1*)EM.DrawBftTDC(1,i+1);
		h[bftnseg+i]->SetAxisRange(	T_min,T_max);
		double PeakPosition = h[bftnseg+i]->GetBinCenter(h[bftnseg+i]->GetMaximumBin());
		double PeakHeight = h[bftnseg+i]->GetMaximum();
		Line[bftnseg+i]= new TLine(PeakPosition,0,PeakPosition,PeakHeight);
		Line[bftnseg+i]->SetLineWidth(2);
		Line[bftnseg+i]->SetLineColor(kBlue);
		Canv_Dn->cd(i+1);
		Param[0]=PeakPosition;
		EM.WriteParameter(ID,Param);
		h[bftnseg+i]->Draw();
		Line[bftnseg+i]->Draw("SAME");
		cout<<i<<endl;
	}
}
void BftSlew(){
	double Tmin = 25,Tmax=60;
	
	TH2* h[2*bftnseg];
	vector<int> ID = {110,0,0,0,2,3};//Cid, PlID, SegID, UorD, Type, nParam For BFT, Set PlID instead of UorD
	vector<double> Param = {0,0,0};
	TF1* func_quad = new TF1("func_quad","[0]*x*x+[1]*x+[2]",Tmin,Tmax);
	func_quad->SetParLimits(0,0,5);

	EM.MakeParameterFile("./Params/BFT_Slew.txt");
	EM.WriteComment("#########################################################");
	EM.WriteComment("# BFT U                                                 #");
	EM.WriteComment("#########################################################");
	
	ID[1]=0;
	TCanvas* Canv_Up1 = new TCanvas("UP1","UP1",1600,1200);
	TCanvas* Canv_Up2 = new TCanvas("UP2","UP2",1600,1200);
	TCanvas* Canv_Up3 = new TCanvas("UP3","UP3",1600,1200);
	TCanvas* Canv_Up4 = new TCanvas("UP4","UP4",1600,1200);
	Canv_Up1->Divide(8,5);
	Canv_Up2->Divide(8,5);
	Canv_Up3->Divide(8,5);
	Canv_Up4->Divide(8,5);

	for(int i=0;i<bftnseg;++i){
		ID[2]=i;
		func_quad->SetRange(Tmin,Tmax);
		if(i<8){
			func_quad->SetRange(35,Tmax);
		}
		h[i]= (TH2*)EM.DrawBftSlewing(0,i+1);
		if(i<40){
			Canv_Up1->cd(i+1);
		}else if(i<80){
			Canv_Up2->cd(i+1-40);
		}else if(i<120){
			Canv_Up3->cd(i+1-80);
		}else{
			Canv_Up4->cd(i+1-120);
		}
		h[i]->Draw("col");
		h[i]->Fit("func_quad","QR");
		double p0 = func_quad->GetParameter(0);
		double p1 = func_quad->GetParameter(1);
		double p2 = func_quad->GetParameter(2);
		Param = {p0,p1,p2};
		EM.WriteParameter(ID,Param);
	}
	
	EM.WriteComment("#########################################################");
	EM.WriteComment("# BFT D                                                 #");
	EM.WriteComment("#########################################################");

	ID[1]=1;
	TCanvas* Canv_Dn1 = new TCanvas("Dn1","Dn1",1600,1200);
	TCanvas* Canv_Dn2 = new TCanvas("Dn2","Dn2",1600,1200);
	TCanvas* Canv_Dn3 = new TCanvas("Dn3","Dn3",1600,1200);
	TCanvas* Canv_Dn4 = new TCanvas("Dn4","Dn4",1600,1200);
	Canv_Dn1->Divide(8,5);
	Canv_Dn2->Divide(8,5);
	Canv_Dn3->Divide(8,5);
	Canv_Dn4->Divide(8,5);
	
	for(int i=0;i<bftnseg;++i){
		ID[2]=i;
		func_quad->SetRange(Tmin,Tmax);
//		if(i<8){
//			func_quad->SetRange(35,Tmax);
//		}
		h[bftnseg+i]= (TH2*)EM.DrawBftSlewing(1,i+1);
		if(i<40){
			Canv_Dn1->cd(i+1);
		}else if(i<80){
			Canv_Dn2->cd(i+1-40);
		}else if(i<120){
			Canv_Dn3->cd(i+1-80);
		}else{
			Canv_Dn4->cd(i+1-120);
		}
		h[bftnseg+i]->Draw("col");
		h[bftnseg+i]->Fit("func_quad","QR");
		double p0 = func_quad->GetParameter(0);
		double p1 = func_quad->GetParameter(1);
		double p2 = func_quad->GetParameter(2);
		Param = {p0,p1,p2};
		EM.WriteParameter(ID,Param);
	}
}



void SchT0(){
	double T_min=450,T_max=500;
	
	TH1* h[schnseg];
	TLine* Line[schnseg];	
	vector<int> ID = {6,0,0,1,0};//Cid, PlID, SegID, AorT, UorD. For BFT, Set PlID instead of UorD
	vector<double> Param = {0,-1};

	EM.MakeParameterFile("./Params/SCH_T0.txt");
	EM.WriteComment("#########################################################");
	EM.WriteComment("# SCH Tdc                                               #");
	EM.WriteComment("#########################################################");
	ID[1]=0;
	TCanvas* Canv_Up = new TCanvas("UP","UP",1600,1200);
	Canv_Up->Divide(8,8);

	for(int i=0;i<schnseg;i++){
		ID[2]=i;
		h[i]= (TH1*)EM.DrawSchTDC(i+1);
		h[i]->SetAxisRange(	T_min,T_max);
		double PeakPosition = h[i]->GetBinCenter(h[i]->GetMaximumBin());
		double PeakHeight = h[i]->GetMaximum();
		Line[i]= new TLine(PeakPosition,0,PeakPosition,PeakHeight);
		Line[i]->SetLineWidth(2);
		Line[i]->SetLineColor(kBlue);
		Canv_Up->cd(i+1);
		Param[0]=PeakPosition;
		EM.WriteParameter(ID,Param);
		h[i]->Draw();
		Line[i]->Draw("SAME");
	}
}
void SchSlew(){
	double Tmin = 40,Tmax=70;
	
	TH2* h[schnseg];
	vector<int> ID = {6,0,0,0,2,3};//Cid, PlID, SegID, UorD, Type, nParam For BFT, Set PlID instead of UorD
	vector<double> Param = {0,0,0};
	TF1* func_quad = new TF1("func_quad","[0]*x*x+[1]*x+[2]",Tmin,Tmax);
	func_quad->SetParLimits(0,0,5);

	EM.MakeParameterFile("./Params/SCH_Slew.txt");
	EM.WriteComment("#########################################################");
	EM.WriteComment("# SCH                                                   #");
	EM.WriteComment("#########################################################");
	
	ID[1]=0;
	TCanvas* Canv_Up1 = new TCanvas("UP1","UP1",1600,1200);
	Canv_Up1->Divide(8,8);

	for(int i=0;i<schnseg;++i){
		ID[2]=i;
		func_quad->SetRange(Tmin,Tmax);
		if(i==63){
			func_quad->SetRange(40,65);
		}
		h[i]= (TH2*)EM.DrawSchSlewing(i+1);
		Canv_Up1->cd(i+1);
		h[i]->Draw("col");
		h[i]->Fit("func_quad","QR");
		double p0 = func_quad->GetParameter(0);
		double p1 = func_quad->GetParameter(1);
		double p2 = func_quad->GetParameter(2);
		Param = {p0,p1,p2};
		EM.WriteParameter(ID,Param);
	}
}
