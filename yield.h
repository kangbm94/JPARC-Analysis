#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "kkana_branch.h"

//#define trigA //all trigA

class yield{

  //#if 0
#if 0
  int condition = 1;
  double vtx[3] = {9.0, 0, -65};
  double vtx_cut[3] = {25,20,100};
  double closeDist_cut = 40;
  double mass_cut[2] = {0.15,0.35};
  //#elif 1
#elif 0
  int condition = 2;
  double vtx[3] = {9.0, 0, -65};
  double vtx_cut[3] = {15,20,100};
  double closeDist_cut = 20;
  double mass_cut[2] = {0.15,0.35};
#else
  int condition = 3;
  double vtx[3] = {9.0, 0, -50};
  double vtx_cut[3] = {30,25,100};
  double closeDist_cut = 100;
  double mass_cut[2] = {0.15,0.40};
#endif

 public :

  class kkana_branch event;

  void sort(TTree *tree, TString name="sort");
  void draw_vertex(TTree *tree, TString name="sort");
  void draw_chisqr(TTree *tree, TString name="sort");
  void kp_yield(TTree *tree, TString name="sort");
  void xi_yield(TTree *tree, TString name="sort");
  void calculate_yield(double p0, double err_p0, double p1, double err_p1,
		       double &yield, double &yield_err);

};

void yield::calculate_yield(double p0, double err_p0, double p1, double err_p1,
			    double &yield, double &yield_err){

  yield = TMath::Sqrt(2*TMath::Pi())*p0*p1;
  yield_err = TMath::Sqrt(2*TMath::Pi())*TMath::Sqrt(p0*p0*err_p1*err_p1+p1*p1*err_p0*err_p0);

  std::cout<<std::endl;
  std::cout<<fixed;
  std::cout.precision(1);
  std::cout<<yield<<" +/- "<<yield_err<<" counts"<<std::endl;

}

void yield::sort(TTree *tree, TString name){

  std::cout<<"condition : "<<condition<<std::endl;
  event.Init(tree);

 //common
 TH1D *btof = new TH1D("btof",";time-of-flight [ns];Counts [/5 ps] ",1000, -2.5, 2.5);

 //chisqr
 TH1D *chisqr_kurama = new TH1D("chisqr_kurama",";ChisqrKurama;Counts",1000, 0, 200);
 TH1D *chisqr_kurama_cut = new TH1D("chisqr_kurama_cut",";ChisqrKurama;Counts",1000, 0, 200);
 TH1D *chisqr_k18 = new TH1D("chisqr_k18",";ChisqrK18;Counts",1000, 0, 30);
 TH1D *chisqr_k18_cut = new TH1D("chisqr_k18_cut",";ChisqrK18;Counts",1000, 0, 10);

 //chisqr cut
 TH1D *vertex_x = new TH1D("vertex_x",";Vertex X [mm];Counts [/1 mm] ",200, -100, 100);
 TH1D *vertex_y = new TH1D("vertex_y",";Vertex Y [mm];Counts [/1 mm] ",200, -100, 100);
 TH1D *vertex_z = new TH1D("vertex_z",";Vertex Z [mm];Counts [/10 mm] ",200, -1000, 1000);
 TH1D *vertex_closedist = new TH1D("vertex_closeDist",";closeDist [mm];Counts [/1 mm] ",100, 0, 100);

 //vertex cut
 TH1D *p_k18 = new TH1D("p_k18",";Momentum [GeV/#font[12]{c}];Counts [/1 MeV/#font[12]{c}] ",400, 1.6, 2.);
 TH1D *p_kurama = new TH1D("p_kurama",";Momentum [GeV/#font[12]{c}];Counts [/20 MeV/#font[12]{c}] ",125, 0, 2.5);
 TH1D *p_kurama_k = new TH1D("p_kurama_k",";Momentum [GeV/#font[12]{c}];Counts [/20 MeV/#font[12]{c}] ",125, 0, 2.5);
 TH1D *vertex_x_cut = new TH1D("vertex_x_cut",";Vertex X [mm];Counts [/1 mm] ",200, -100, 100);
 TH1D *vertex_y_cut = new TH1D("vertex_y_cut",";Vertex Y [mm];Counts [/1 mm] ",200, -100, 100);
 TH1D *vertex_z_cut = new TH1D("vertex_z_cut",";Vertex Z [mm];Counts [/10 mm] ",200, -1000, 1000);
 TH1D *vertex_closedist_cut = new TH1D("vertex_closeDist_cut",";closeDist [mm];Counts [/1 mm] ",100, 0, 100);
 TH2D *p_m = new TH2D("p_m",";Mass/Charge [GeV/#font[12]{c}^{2}/q];Momentum [GeV/#font[12]{c}]",400,-1,2,400,0,2.5);

 //momentum cut
 TH1D *mass = new TH1D("mass",";Mass/Charge [GeV/#font[12]{c}^{2}/q];Counts [/10 MeV/#font[12]{c}^{2}/q] ",280, -0.8, 2);

 //k+ cut
 TH1D *mm = new TH1D("mm",";Missing Mass [GeV/#font[12]{c}^{2}];Counts [/5 MeV/#font[12]{c}^{2}] ",160, 1., 1.8);
 TH1D *m2_kp = new TH1D("m2_kp",";Mass^{2}/Charge [(GeV/#font[12]{c}^{2})^{2}/q];Counts [/0.001 (GeV/#font[12]{c}^{2})^{2}/q ",400, 0.1, 0.5);
 TH1D *p_kp = new TH1D("p_kp",";Momentum [GeV/#font[12]{c}];Counts [/5 MeV/#font[12]{c}] ",400, 0, 2.);
 TH1D *tracks = new TH1D("tracks",";Tracks;Counts", 10, 0, 10);

 TFile* newtf = new TFile(Form("%s.root",name.Data()),"recreate");

 int entry = tree -> GetEntries();
 std::cout<<"# of events : "<<entry<<std::endl;
 for(int i=0;i<entry;i++){
   tree -> GetEntry(i);
   if(i%10000==0) std::cout<<i<<" / "<<entry<<std::endl;

   //common
   btof -> Fill(event.CBtof0);

#ifdef trigA
   //trig.A
   if(event.trigflag[20]<0 || event.trigpat[20]!=21) continue;
#endif
   int ikk = -1; int ikk_cut = 0; int ibeam=0; int flag=-1;
   for(int ikm=0;ikm<event.nKm;ikm++){
     for(int ikp=0;ikp<event.nKp;ikp++){
       ikk++;

       if(event.nKK==1&&ikm==0) chisqr_k18 -> Fill(event.chisqrK18[ikm]);
       if(event.nKK==1&&ikp==0) chisqr_kurama -> Fill(event.chisqrKurama[ikp]);

       //chisqr cut
       if(event.chisqrK18[ikm] > 20 || event.chisqrKurama[ikp] > 200) continue;

       if(event.nKK==1){
	 chisqr_k18_cut -> Fill(event.chisqrK18[ikm]);
	 chisqr_kurama_cut -> Fill(event.chisqrKurama[ikp]);
       }
       vertex_x -> Fill(event.vtx[ikk]);
       vertex_y -> Fill(event.vty[ikk]);
       vertex_z -> Fill(event.vtz[ikk]);
       vertex_closedist -> Fill(event.closeDist[ikk]);

       //vertex cut
       int vtx_flag = 0;
       if(TMath::Abs(event.vtx[ikk]-vtx[0])>vtx_cut[0]) vtx_flag++;
       if(TMath::Abs(event.vty[ikk]-vtx[1])>vtx_cut[1]) vtx_flag++;
       if(TMath::Abs(event.vtz[ikk]-vtx[2])>vtx_cut[2]) vtx_flag++;
       if(event.closeDist[ikk]>closeDist_cut) vtx_flag++;
       if(vtx_flag!=0) continue;

       vertex_x_cut -> Fill(event.vtx[ikk]);
       vertex_y_cut -> Fill(event.vty[ikk]);
       vertex_z_cut -> Fill(event.vtz[ikk]);
       vertex_closedist_cut -> Fill(event.closeDist[ikk]);

       p_k18 -> Fill(event.pK18[ikm]);
       p_kurama -> Fill(event.pKurama[ikp]);
       p_m -> Fill(event.qKurama[ikp]*TMath::Sqrt(event.m2[ikp]),event.pKurama[ikp]);


       ikk_cut++;

       //momentum cut
       if(event.pKurama[ikp]>1.4) continue;

       mass -> Fill(TMath::Sqrt(event.m2[ikp])*event.qKurama[ikp]);

       //k+ cut
       if(event.qKurama[ikp]<0) continue;
       if(event.m2[ikp]<mass_cut[0]||event.m2[ikp]>mass_cut[1]) continue;

       if(flag!=ikm){
	 ibeam++;
	 flag=ikm;
       }

       mm -> Fill(event.MissMassCorr[ikp]);
       m2_kp -> Fill(event.m2[ikp]*event.qKurama[ikp]);
       p_kp -> Fill(event.pKurama[ikp]);
       p_kurama_k -> Fill(event.pKurama[ikp]);

     }
     if(ibeam==1&&ikk_cut!=0) tracks -> Fill(ikk_cut);
   }
 }
 std::cout<<"tracks : "<<tracks -> GetEntries()<<std::endl;

 newtf -> cd();
 btof -> Write();
 chisqr_k18 -> Write();
 chisqr_kurama -> Write();

 //chisqr cut
 chisqr_k18_cut -> Write();
 chisqr_kurama_cut -> Write();
 vertex_x -> Write();
 vertex_y -> Write();
 vertex_z -> Write();
 vertex_closedist -> Write();

 //vertex cut
 p_k18 -> Write();
 p_kurama -> Write();
 p_kurama_k -> Write();
 vertex_x_cut -> Write();
 vertex_y_cut -> Write();
 vertex_z_cut -> Write();
 vertex_closedist_cut -> Write();
 p_m -> Write();

 //momentum cut
 mass -> Write();
 mm -> Write();
 p_kp -> Write();
 m2_kp -> Write();
 tracks -> Write();

 newtf -> Close();
 delete newtf;

}

void yield::draw_vertex(TTree *tree, TString name){

  TFile* tf = new TFile(Form("%s.root",name.Data()),"read");
  TCanvas *can_vertex = new TCanvas("can_vertex","can_vertex",1200,600);
  can_vertex -> Divide(2,2);
  TCanvas *can_track = new TCanvas("can_track","can_track",1200,600);

  TH1D *vertex_x_cut = (TH1D *) tf -> Get("vertex_x_cut");
  TH1D *vertex_y_cut = (TH1D *) tf -> Get("vertex_y_cut");
  TH1D *vertex_z_cut = (TH1D *) tf -> Get("vertex_z_cut");
  TH1D *vertex_closeDist_cut = (TH1D *) tf -> Get("vertex_closeDist_cut");

  TH1D *vertex_x = (TH1D *) tf -> Get("vertex_x");
  TH1D *vertex_y = (TH1D *) tf -> Get("vertex_y");
  TH1D *vertex_z = (TH1D *) tf -> Get("vertex_z");
  TH1D *vertex_closeDist = (TH1D *) tf -> Get("vertex_closeDist");

  TH1D *tracks = (TH1D *) tf -> Get("tracks");

  can_vertex -> cd(1);
  vertex_x -> Draw();
  vertex_x_cut -> Draw("same");
  vertex_x_cut -> SetLineColor(2);
  can_vertex -> cd(2);
  vertex_y -> Draw();
  vertex_y_cut -> Draw("same");
  vertex_y_cut -> SetLineColor(2);
  can_vertex -> cd(3);
  vertex_z -> Draw();
  vertex_z_cut -> Draw("same");
  vertex_z_cut -> SetLineColor(2);
  can_vertex -> cd(4);
  vertex_closeDist -> Draw();
  vertex_closeDist_cut -> Draw("same");
  vertex_closeDist_cut -> SetLineColor(2);

  can_track -> cd();
  tracks -> Draw();

  can_vertex -> SaveAs(Form("vertex_%s.pdf",name.Data()));
  can_track -> SaveAs(Form("track_%s.pdf",name.Data()));

}

void yield::draw_chisqr(TTree *tree, TString name){

  TFile* tf = new TFile(Form("%s.root",name.Data()),"read");
  TCanvas *can_chisqr = new TCanvas("can_chisqr","can_chisqr",1200,600);
  can_chisqr -> Divide(2,2);

  TH1D *chisqr_kurama = (TH1D *) tf -> Get("chisqr_kurama");
  TH1D *chisqr_kurama_cut = (TH1D *) tf -> Get("chisqr_kurama_cut");
  TH1D *chisqr_k18 = (TH1D *) tf -> Get("chisqr_k18");
  TH1D *chisqr_k18_cut = (TH1D *) tf -> Get("chisqr_k18_cut");

  TH1D *chisqr_kurama2 = (TH1D *) chisqr_kurama -> Clone();
  TH1D *chisqr_kurama_cut2 = (TH1D *) chisqr_kurama_cut -> Clone();
  TH1D *chisqr_k182 = (TH1D *) chisqr_k18 -> Clone();
  TH1D *chisqr_k18_cut2 = (TH1D *) chisqr_k18_cut -> Clone();

  can_chisqr -> cd(1);
  //gPad -> SetLogy();
  chisqr_kurama -> Draw();
  chisqr_kurama_cut -> Draw("same");
  chisqr_kurama_cut -> SetLineColor(2);

  can_chisqr -> cd(2);
  //gPad -> SetLogy();
  chisqr_kurama2 -> Draw();
  chisqr_kurama2 -> GetXaxis() -> SetRangeUser(0,20);
  chisqr_kurama_cut2 -> Draw("same");
  chisqr_kurama_cut2 -> SetLineColor(2);
  chisqr_kurama_cut2 -> GetXaxis() -> SetRangeUser(0,20);

  can_chisqr -> cd(3);
  //gPad -> SetLogy();
  chisqr_k18 -> Draw();
  chisqr_k18_cut -> Draw("same");
  chisqr_k18_cut -> SetLineColor(2);

  can_chisqr -> cd(4);
  //gPad -> SetLogy();
  chisqr_k182 -> Draw();
  chisqr_k182 -> GetXaxis() -> SetRangeUser(0,5);
  chisqr_k18_cut2 -> Draw("same");
  chisqr_k18_cut2 -> SetLineColor(2);
  chisqr_k18_cut2 -> GetXaxis() -> SetRangeUser(0,5);

  can_chisqr -> SaveAs(Form("chisqr_%s.pdf",name.Data()));

}

void yield::kp_yield(TTree *tree, TString name){

  TFile* tf = new TFile(Form("%s.root",name.Data()),"read");
  TCanvas *can_kp = new TCanvas("can_kp","can_kp",1200,600);

  double par[5];
  TF1 *func_yield = new TF1("func_yield","gaus(0)+[3]*x+[4]",0.15,0.35);
  TF1 *func_kp = new TF1("func_kp","gaus(0)",0.15,0.35);
  TF1 *func_bkg = new TF1("func_bkg","[0]*x+[1]",0.15,0.35);
  TH1D *m2_kp = (TH1D *) tf -> Get("m2_kp");
  can_kp -> cd();
  m2_kp -> Draw();
  m2_kp -> Fit("func_yield","REQ");
  func_yield -> GetParameters(par);
  func_kp -> SetParameters(&par[0]);
  func_bkg -> SetParameters(&par[3]);
  func_kp -> Draw("same");
  func_kp -> SetLineColor(3);
  func_bkg -> Draw("same");
  func_bkg -> SetLineColor(1);
  func_bkg -> SetLineStyle(2);

  can_kp -> SaveAs(Form("kp_yield_%s.pdf",name.Data()));

}

void yield::xi_yield(TTree *tree, TString name){

  TFile* tf = new TFile(Form("%s.root",name.Data()),"read");
  TCanvas *can_xi = new TCanvas("can_xi","can_xi",1200,600);

  double par[6]={0}; double err_par[6]={0};
  TF1 *func_yield = new TF1("func_yield","gaus(0)+gaus(3)",1.28,1.4);
  TF1 *func_xi = new TF1("func_xi","gaus(0)",1.28,1.4);
  TF1 *func_bkg = new TF1("func_bkg","gaus(0)",1.28,1.4);

  TH1D *mm = (TH1D *) tf -> Get("mm");
  can_xi -> cd();
  mm -> Draw();
  double height = mm -> GetMaximum();
  double xi_center = mm -> GetBinCenter(mm -> GetMaximumBin());
  func_yield -> SetParLimits(0,height*0.5,height);
  func_yield -> SetParLimits(1,xi_center-0.3,xi_center+0.3);
  func_yield -> SetParLimits(2,0.01,0.03);
  func_yield -> SetParLimits(3,height*0.,height*0.5);
  func_yield -> SetParLimits(4,xi_center-0.3,xi_center+0.3);
  func_yield -> SetParLimits(5,0.02,0.10);
  mm -> Fit("func_yield","RE");
  func_yield -> GetParameters(par);
  err_par[0] = func_yield -> GetParError(0);
  err_par[1] = func_yield -> GetParError(1);

  func_xi -> SetParameters(&par[0]);
  func_bkg -> SetParameters(&par[3]);
  func_xi -> Draw("same");
  func_xi -> SetLineColor(3);
  func_bkg -> Draw("same");
  func_bkg -> SetLineColor(1);
  func_bkg -> SetLineStyle(2);

  can_xi -> SaveAs(Form("xi_yield_%s.pdf",name.Data()));
  double xi_yield; double err_xi_yield;
  std::cout<<"total entries : "<<  mm -> GetEntries()<<std::endl;
  calculate_yield(par[0], err_par[0], par[1], err_par[1], xi_yield, err_xi_yield);

}
