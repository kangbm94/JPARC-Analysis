double dgaus(double* x, double* p){
	double val=0;
	val=Gaussian(x[0],p[0],p[1],p[2])+Gaussian(x[0],p[0]+p[3],p[4],p[5]);//Mean,sig, amp
	return val;
}
TF1* func_dgaus=new TF1("func_dgaus","dgaus",400,1000,6);
void DgausFit(){
	TString fn = "HistMTSdcInParameter05340.root";
	TFile* file = new TFile(fn,"read");
	TString hnbase = "MeanTime";
	TH1D* hist[5];
	double T0RangeMin[10]={430,430,430,430,430,430,930,930,940,940};
	double T0RangeMax[10]={450,450,450,450,450,450,960,960,960,960};
	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(3,2);
	double cut[5];
	for(int i=0;i<5;i++){
		c1->cd(i+1);
		int l1 = 2*i,l2=2*i+1;
		double r1=T0RangeMin[l1]-40,r2=T0RangeMin[l1]+10;
		func_dgaus->SetParLimits(0,r1-10,T0RangeMax[l1]);
		func_dgaus->SetParLimits(1,2,25);
		func_dgaus->SetParLimits(3,-30,-10);
		func_dgaus->SetParLimits(4,2,20);
		TString hn =hnbase+Form("L%dL%d",l1,l2);
		hist[i]=(TH1D*)file->Get(hn);
		hist[i]->Fit("func_dgaus","QM");
		double mean = func_dgaus->GetParameter(0);
		double dt = func_dgaus->GetParameter(3);
		cut[i]= mean+dt/2;
		cout<<"Layer "<<l1<<" : "<<cut[i]<<endl;
	}
}
