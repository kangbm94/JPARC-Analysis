#include "ToFMethod.hh"
int Run = 5426;
vector<int> Exception={1,2,3,5,7,8,9,12,14};

vector<double> Origin = {4,0,-67};
vector<double> TargetSize = {25,20,65};



void ToF(){
	int runnum = 5641;
	t.LoadFile("./rootfiles/CH2/AllKuramaHodoscope.root");
//	t.LoadFile(Form("./rootfiles/CH2/DstKuramaHodoscope0%d.root",runnum));
//	t.LoadFile("./rootfiles/Dsts/DstKuramaHodoscope05395.root");
//	t.LoadFile("./rootfiles/Pion/DstKuramaHodoscope05323.root");
//	t.LoadFile("~/WS_data/ch2target/run05641_KuramaHodoscopeTest.root");
	t.LoadKHodo();
}

void FitADC(){
	TH1* PEDUp[tofnseg];
	TH1* PEDDn[tofnseg];
	TH1* ADCUp[tofnseg];
	TH1* ADCDn[tofnseg];
	TCanvas* Can_PedU = new TCanvas("Can_PedU","Can_PedU",1500,700);
	Can_PedU->Divide(6,4);
	TCanvas* Can_PedD = new TCanvas("Can_PedD","Can_PedD",1500,700);
	Can_PedD->Divide(6,4);
	TCanvas* Can_AdcU = new TCanvas("Can_AdcU","Can_AdcU",1500,700);
	Can_AdcU->Divide(6,4);
	TCanvas* Can_AdcD = new TCanvas("Can_AdcD","Can_AdcD",1500,700);
	Can_AdcD->Divide(6,4);
	double Base_Up[24];
	double Base_Dn[24];
	double MPV_Up[24];
	double MPV_Dn[24];
	
	for(int i=0;i<tofnseg;++i){
		int seg = i+1;
		PEDUp[i]=t.GetPEDHisto(seg,0);
		PEDDn[i]=t.GetPEDHisto(seg,1);
		ADCUp[i]=t.GetADCHisto(seg,0,0);
		ADCDn[i]=t.GetADCHisto(seg,1,0);
		if(seg>6){
			ADCUp[i]=t.GetADCHisto(seg,0,0);
			ADCDn[i]=t.GetADCHisto(seg,1,0);
		}
		PEDUp[i]->SetAxisRange(0,300);
		PEDDn[i]->SetAxisRange(0,300);
		ADCUp[i]->SetAxisRange(500,1200);
		ADCDn[i]->SetAxisRange(500,1200);
		if(seg==1){ADCDn[i]=t.GetADCHisto(seg,1,0);ADCDn[i]->SetAxisRange(600,1000);}
		if(seg==1){ADCUp[i]=t.GetADCHisto(seg,0,0);ADCUp[i]->SetAxisRange(600,1000);}
		if(seg==2){ADCUp[i]=t.GetADCHisto(seg,0,1);ADCUp[i]->SetAxisRange(500,1000);}
		if(seg==2){ADCDn[i]=t.GetADCHisto(seg,1,1);ADCDn[i]->SetAxisRange(400,1000);}
		if(seg==4){ADCUp[i]->SetAxisRange(150,250);}
		if(seg==9){PEDDn[i]->SetAxisRange(0,100);}
		Can_PedU->cd(i+1);
		PEDUp[i]->Fit("f_gaus");
		Base_Up[i]=f_gaus->GetParameter(1);
		Can_PedD->cd(i+1);
		PEDDn[i]->Fit("f_gaus");
		Base_Dn[i]=f_gaus->GetParameter(1);
		Can_AdcU->cd(i+1);
		ADCUp[i]->Fit("f_landau");
		MPV_Up[i]=f_landau->GetParameter(1);
		Can_AdcD->cd(i+1);
		ADCDn[i]->Fit("f_landau");
		MPV_Dn[i]=f_landau->GetParameter(1);
	}
	t.MakeParameterFile("./param/ToF_ADC.txt");
	
	vector<int>ID={7,0,0,0,0};
	vector<double>Param={0,0};
	
	t.WriteComment("#########################################################");
	t.WriteComment(Form("# TOF Adc Up, Callibrated with run 0%d                #",Run));
	t.WriteComment("#########################################################");

	for(int i=0;i<tofnseg;++i){
		ID[2]=i;
		Param[0]=Base_Up[i];
		Param[1]=MPV_Up[i];
		t.WriteParameter(ID,Param);
	}
	ID[4]=1;
	t.WriteComment("#########################################################");
	t.WriteComment("# TOF Adc Down                                          #");
	t.WriteComment("#########################################################");
	for(int i=0;i<tofnseg;++i){
		ID[2]=i;
		Param[0]=Base_Dn[i];
		Param[1]=MPV_Dn[i];
		t.WriteParameter(ID,Param);
	}
}
void ViewPHC(int seg,int particle=0){
	TCanvas* c1 = new TCanvas("c1","c1",1500,700);
	c1->Divide(2,2);
	TH2* h[4];
	h[0]=t.GetPHCHisto(seg,0,particle);
	h[1]=t.GetPHCedHisto(seg,0,particle);
	h[2]=t.GetPHCHisto(seg,1,particle);
	h[3]=t.GetPHCedHisto(seg,1,particle);
	for(int i=0;i<4;i++){
		h[i]->SetOption("col");
		c1->cd(i+1);
		gPad->SetLogz();
		h[i]->Draw();
	}
}


void RestoreCallib(){
	ToFManager Org;
	ToFManager Wrong;
	ToFManager Correct;

	double OrgVTOF[tofnseg];
	double WrongVTOF[tofnseg];
	double CorrectVTOF[tofnseg];
	Org.LoadOldVTOF("VTOF_before",OrgVTOF);
	Wrong.LoadOldVTOF("VTOF_after",WrongVTOF);
	Correct.MakeParameterFile("VTOF_correct");
	vector<int>ID = {7,0,0,1,2};
	vector<double>Param={0,-1};
	for(int i=0; i<tofnseg; ++i){
		ID[2]=i;
		Param[0]=2*OrgVTOF[i]-WrongVTOF[i];
		Correct.WriteParameter(ID,Param);
	}

}


void CallibratePion(int conf=0){
	if(conf){
		CallibrateParticleFromChain(1);
	}
	else{
		CallibrateParticleFromHist(1);
	}
}
void CallibrateKaon(int conf=0){
	if(conf){
		CallibrateParticleFromChain(2);
	}
	else{
		CallibrateParticleFromHist(2);
	}
}
void CallibrateProton(int conf = 0){
	if(conf){
		CallibrateParticleFromChain(3);
	}
	else{
		CallibrateParticleFromHist(3);
	}
}


void DoPHC(){
	gStyle->SetOptFit(1111);
	t.MakeParameterFile("./param/ToF_PHC.txt");
	int particle = 2;//0->Pion,1->Kaon, 2->Proton
	int mod = 1;
	TCanvas* Canv_Up = new TCanvas("Canv_Up","Canv_Up",1500,700);
	Canv_Up->Divide(6,4);
	TCanvas* Canv_Dn = new TCanvas("Canv_Dn","Canv_Dn",1500,700);
	Canv_Dn->Divide(6,4);
	double p0_up[24];double p1_up[24];double p2_up[24];
	double p0_dn[24];double p1_dn[24];double p2_dn[24];

	vector<int>ID={7,0,0,0,0,3};
	vector<double>Param={0,0,0};
	
	t.WriteComment("#########################################################");
	t.WriteComment(Form("# TOF PHC Up, Callibrated with run 0%d                #",Run));
	t.WriteComment("#########################################################");
	mod=0;

	for(int i=0;i<tofnseg;++i){
		int seg = i+1;
		ID[2]=i;
		if(seg==1){
			ID[4]=0;
		}else{
			ID[4]=1;
		}
		Canv_Up->cd(i+1);
		if(seg<7){
			particle =0;
		}
		else{
			particle=0;
		}
		slew_func->SetRange(0,5);
		if(seg<20){
//			mod=1;
		}
		else{
		//	mod=0;
		}
		t.FitTimewalk(seg,0,particle,mod);
		p0_up[i]=slew_func->GetParameter(0);
		p1_up[i]=slew_func->GetParameter(1);
		p2_up[i]=slew_func->GetParameter(2);
//		p2_up[i]+=SlewFunc(1,p0_up[i],p1_up[i],p2_up[i]);
//		cout<<SlewFunc(1,p0_up[i],p1_up[i],p2_up[i])<<endl;
		Param[0]=-p0_up[i];
		Param[1]=p1_up[i];
		Param[2]=-p2_up[i];
		t.WriteParameter(ID,Param);
	}
	
	t.WriteComment("#########################################################");
	t.WriteComment(Form("# TOF PHC Dn, Callibrated with run 0%d                #",Run));
	t.WriteComment("#########################################################");
	ID[3]=1;	
	for(int i=0;i<tofnseg;++i){
		int seg = i+1;
		ID[2]=i;
		if(seg==1){
			ID[4]=0;
		}else{
			ID[4]=1;
		}
		Canv_Dn->cd(i+1);
		if(seg<7){
			particle =0;
		}
		else{
			particle=0;
		}
		if(seg<20){
//			mod=1;
		}
		else{
		//	mod=0;
		}
		t.FitTimewalk(seg,1,particle,mod);
		p0_dn[i]=slew_func->GetParameter(0);
		p1_dn[i]=slew_func->GetParameter(1);
		p2_dn[i]=slew_func->GetParameter(2);
//		p2_dn[i]+=SlewFunc(1,p0_dn[i],p1_dn[i],p2_dn[i]);
		
		Param[0]=-p0_dn[i];
		Param[1]=p1_dn[i];
		Param[2]=-p2_dn[i];
		t.WriteParameter(ID,Param);
	}
}
void DoPHC(int seg,int UD){
	gStyle->SetOptFit(1111);
	t.SaveHisto(Form("ToFPHC_%d.root",seg));
	t.FitTimewalk(seg,UD,-1,0);
}

void FitTDC(){
}
