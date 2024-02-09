//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 26 18:10:32 2023 by ROOT version 6.26/02
// from TTree kk/tree of KkAna
// found on file: DstHSKKAna.root
//////////////////////////////////////////////////////////

#ifndef HSKKAna_h
#define HSKKAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class HSKKAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TTree          *TreeOut; 
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runnum;
   Int_t           evnum;
   Int_t           spill;
   Int_t           trigpat[32];
   Int_t           trigflag[32];
   Int_t           nhBh1;
   Int_t           csBh1[5];   //[nhBh1]
   Double_t        Bh1Seg[5];   //[nhBh1]
   Double_t        tBh1[5];   //[nhBh1]
   Double_t        dtBh1[5];   //[nhBh1]
   Double_t        deBh1[5];   //[nhBh1]
   Int_t           nhBh2;
   Int_t           csBh2[11];   //[nhBh2]
   Double_t        Bh2Seg[11];   //[nhBh2]
   Double_t        tBh2[11];   //[nhBh2]
   Double_t        t0Bh2[11];   //[nhBh2]
   Double_t        dtBh2[11];   //[nhBh2]
   Double_t        deBh2[11];   //[nhBh2]
   Double_t        Time0Seg;
   Double_t        deTime0;
   Double_t        Time0;
   Double_t        CTime0;
   Double_t        Btof0Seg;
   Double_t        deBtof0;
   Double_t        Btof0;
   Double_t        CBtof0;
   Int_t           nhTof;
   Int_t           csTof[26];   //[nhTof]
   Double_t        TofSeg[26];   //[nhTof]
   Double_t        tTof[26];   //[nhTof]
   Double_t        dtTof[26];   //[nhTof]
   Double_t        deTof[26];   //[nhTof]
   Int_t           nhBvh;
   Double_t        BvhSeg[8];   //[nhBvh]
   Int_t           nhBft;
   Int_t           csBft[19];   //[nhBft]
   Double_t        tBft[19];   //[nhBft]
   Double_t        wBft[19];   //[nhBft]
   Double_t        BftPos[19];   //[nhBft]
   Double_t        BftSeg[19];   //[nhBft]
   Int_t           nhSch;
   Int_t           csSch[36];   //[nhSch]
   Double_t        tSch[36];   //[nhSch]
   Double_t        wSch[36];   //[nhSch]
   Double_t        SchPos[36];   //[nhSch]
   Double_t        SchSeg[36];   //[nhSch]
   Int_t           nlBcOut;
   Int_t           ntBcOut;
   Int_t           nhBcOut[6];   //[ntBcOut]
   Double_t        chisqrBcOut[6];   //[ntBcOut]
   Double_t        x0BcOut[6];   //[ntBcOut]
   Double_t        y0BcOut[6];   //[ntBcOut]
   Double_t        u0BcOut[6];   //[ntBcOut]
   Double_t        v0BcOut[6];   //[ntBcOut]
   Double_t        xtgtBcOut[6];   //[ntBcOut]
   Double_t        ytgtBcOut[6];   //[ntBcOut]
   Double_t        xbh2BcOut[6];   //[ntBcOut]
   Double_t        ybh2BcOut[6];   //[ntBcOut]
   Int_t           ntK18;
   Int_t           nhK18[3];   //[ntK18]
   Double_t        chisqrK18[3];   //[ntK18]
   Double_t        pK18[3];   //[ntK18]
   Double_t        xtgtHS[3];   //[ntK18]
   Double_t        ytgtHS[3];   //[ntK18]
   Double_t        utgtHS[3];   //[ntK18]
   Double_t        vtgtHS[3];   //[ntK18]
   Double_t        thetaHS[3];   //[ntK18]
   Double_t        tpcHSvpx[3][4];   //[ntK18]
   Double_t        tpcHSvpy[3][4];   //[ntK18]
   Double_t        tpcHSvpz[3][4];   //[ntK18]
   Int_t           nlSdcIn;
   Int_t           ntSdcIn;
   Int_t           nhSdcIn[8];   //[ntSdcIn]
   Double_t        chisqrSdcIn[8];   //[ntSdcIn]
   Double_t        x0SdcIn[8];   //[ntSdcIn]
   Double_t        y0SdcIn[8];   //[ntSdcIn]
   Double_t        u0SdcIn[8];   //[ntSdcIn]
   Double_t        v0SdcIn[8];   //[ntSdcIn]
   Int_t           nlSdcOut;
   Int_t           ntSdcOut;
   Int_t           nhSdcOut[8];   //[ntSdcOut]
   Double_t        chisqrSdcOut[8];   //[ntSdcOut]
   Double_t        x0SdcOut[8];   //[ntSdcOut]
   Double_t        y0SdcOut[8];   //[ntSdcOut]
   Double_t        u0SdcOut[8];   //[ntSdcOut]
   Double_t        v0SdcOut[8];   //[ntSdcOut]
   Int_t           ntKurama;
   Int_t           nhKurama[15];   //[ntKurama]
   Double_t        chisqrKurama[15];   //[ntKurama]
   Double_t        stof[15];   //[ntKurama]
   Double_t        cstof[15];   //[ntKurama]
   Double_t        path[15];   //[ntKurama]
   Double_t        pKurama[15];   //[ntKurama]
   Double_t        qKurama[15];   //[ntKurama]
   Double_t        m2[15];   //[ntKurama]
   Double_t        xtgtKurama[15];   //[ntKurama]
   Double_t        ytgtKurama[15];   //[ntKurama]
   Double_t        utgtKurama[15];   //[ntKurama]
   Double_t        vtgtKurama[15];   //[ntKurama]
   Double_t        thetaKurama[15];   //[ntKurama]
   Double_t        xtofKurama[15];   //[ntKurama]
   Double_t        ytofKurama[15];   //[ntKurama]
   Double_t        utofKurama[15];   //[ntKurama]
   Double_t        vtofKurama[15];   //[ntKurama]
   Double_t        tofsegKurama[15];   //[ntKurama]
   Double_t        best_deTof[15];   //[ntKurama]
   Double_t        best_TofSeg[15];   //[ntKurama]
   Double_t        vpxtpc[15][5];   //[ntKurama]
   Double_t        vpytpc[15][5];   //[ntKurama]
   Double_t        vpztpc[15][5];   //[ntKurama]
   Int_t           nKm;
   Int_t           nKp;
   Int_t           nKK;
   Double_t        vtx[24];   //[nKK]
   Double_t        vty[24];   //[nKK]
   Double_t        vtz[24];   //[nKK]
   Double_t        KMPX[24];   //[nKK]
   Double_t        KMPY[24];   //[nKK]
   Double_t        KMPZ[24];   //[nKK]
   Double_t        KPPX[24];   //[nKK]
   Double_t        KPPY[24];   //[nKK]
   Double_t        KPPZ[24];   //[nKK]
   Double_t        closeDist[24];   //[nKK]
   Int_t           inside[24];   //[nKK]
   Double_t        theta[24];   //[nKK]
   Double_t        MissMass[24];   //[nKK]
   Double_t        MissMassCorr[24];   //[nKK]
   Double_t        MissMassCorrDE[24];   //[nKK]
   Double_t        MissMomx[24];   //[nKK]
   Double_t        MissMomy[24];   //[nKK]
   Double_t        MissMomz[24];   //[nKK]
   Double_t        MissMomxCor[24];   //[nKK]
   Double_t        MissMomyCor[24];   //[nKK]
   Double_t        MissMomzCor[24];   //[nKK]
   Double_t        thetaCM[24];   //[nKK]
   Double_t        costCM[24];   //[nKK]
   Int_t           Kflag[24];   //[nKK]
   Bool_t          Xiflag[24];   //[nKK]
   Bool_t          KKflag[24];   //[nKK]
   Double_t        xkm[24];   //[nKK]
   Double_t        ykm[24];   //[nKK]
   Double_t        ukm[24];   //[nKK]
   Double_t        vkm[24];   //[nKK]
   Double_t        xkp[24];   //[nKK]
   Double_t        ykp[24];   //[nKK]
   Double_t        ukp[24];   //[nKK]
   Double_t        vkp[24];   //[nKK]
   Double_t        pOrg[24];   //[nKK]
   Double_t        pCalc[24];   //[nKK]
   Double_t        pCorr[24];   //[nKK]
   Double_t        pCorrDE[24];   //[nKK]

   // List of branches
   TBranch        *b_runnum;   //!
   TBranch        *b_evnum;   //!
   TBranch        *b_spill;   //!
   TBranch        *b_trigpat;   //!
   TBranch        *b_trigflag;   //!
   TBranch        *b_nhBh1;   //!
   TBranch        *b_csBh1;   //!
   TBranch        *b_Bh1Seg;   //!
   TBranch        *b_tBh1;   //!
   TBranch        *b_dtBh1;   //!
   TBranch        *b_deBh1;   //!
   TBranch        *b_nhBh2;   //!
   TBranch        *b_csBh2;   //!
   TBranch        *b_Bh2Seg;   //!
   TBranch        *b_tBh2;   //!
   TBranch        *b_t0Bh2;   //!
   TBranch        *b_dtBh2;   //!
   TBranch        *b_deBh2;   //!
   TBranch        *b_Time0Seg;   //!
   TBranch        *b_deTime0;   //!
   TBranch        *b_Time0;   //!
   TBranch        *b_CTime0;   //!
   TBranch        *b_Btof0Seg;   //!
   TBranch        *b_deBtof0;   //!
   TBranch        *b_Btof0;   //!
   TBranch        *b_CBtof0;   //!
   TBranch        *b_nhTof;   //!
   TBranch        *b_csTof;   //!
   TBranch        *b_TofSeg;   //!
   TBranch        *b_tTof;   //!
   TBranch        *b_dtTof;   //!
   TBranch        *b_deTof;   //!
   TBranch        *b_nhBvh;   //!
   TBranch        *b_BvhSeg;   //!
   TBranch        *b_nhBft;   //!
   TBranch        *b_csBft;   //!
   TBranch        *b_tBft;   //!
   TBranch        *b_wBft;   //!
   TBranch        *b_BftPos;   //!
   TBranch        *b_BftSeg;   //!
   TBranch        *b_nhSch;   //!
   TBranch        *b_csSch;   //!
   TBranch        *b_tSch;   //!
   TBranch        *b_wSch;   //!
   TBranch        *b_SchPos;   //!
   TBranch        *b_SchSeg;   //!
   TBranch        *b_nlBcOut;   //!
   TBranch        *b_ntBcOut;   //!
   TBranch        *b_nhBcOut;   //!
   TBranch        *b_chisqrBcOut;   //!
   TBranch        *b_x0BcOut;   //!
   TBranch        *b_y0BcOut;   //!
   TBranch        *b_u0BcOut;   //!
   TBranch        *b_v0BcOut;   //!
   TBranch        *b_xtgtBcOut;   //!
   TBranch        *b_ytgtBcOut;   //!
   TBranch        *b_xbh2BcOut;   //!
   TBranch        *b_ybh2BcOut;   //!
   TBranch        *b_ntK18;   //!
   TBranch        *b_nhK18;   //!
   TBranch        *b_chisqrK18;   //!
   TBranch        *b_pK18;   //!
   TBranch        *b_xtgtHS;   //!
   TBranch        *b_ytgtHS;   //!
   TBranch        *b_utgtHS;   //!
   TBranch        *b_vtgtHS;   //!
   TBranch        *b_thetaHS;   //!
   TBranch        *b_tpcHSvpx;   //!
   TBranch        *b_tpcHSvpy;   //!
   TBranch        *b_tpcHSvpz;   //!
   TBranch        *b_nlSdcIn;   //!
   TBranch        *b_ntSdcIn;   //!
   TBranch        *b_nhSdcIn;   //!
   TBranch        *b_chisqrSdcIn;   //!
   TBranch        *b_x0SdcIn;   //!
   TBranch        *b_y0SdcIn;   //!
   TBranch        *b_u0SdcIn;   //!
   TBranch        *b_v0SdcIn;   //!
   TBranch        *b_nlSdcOut;   //!
   TBranch        *b_ntSdcOut;   //!
   TBranch        *b_nhSdcOut;   //!
   TBranch        *b_chisqrSdcOut;   //!
   TBranch        *b_x0SdcOut;   //!
   TBranch        *b_y0SdcOut;   //!
   TBranch        *b_u0SdcOut;   //!
   TBranch        *b_v0SdcOut;   //!
   TBranch        *b_ntKurama;   //!
   TBranch        *b_nhKurama;   //!
   TBranch        *b_chisqrKurama;   //!
   TBranch        *b_stof;   //!
   TBranch        *b_cstof;   //!
   TBranch        *b_path;   //!
   TBranch        *b_pKurama;   //!
   TBranch        *b_qKurama;   //!
   TBranch        *b_m2;   //!
   TBranch        *b_xtgtKurama;   //!
   TBranch        *b_ytgtKurama;   //!
   TBranch        *b_utgtKurama;   //!
   TBranch        *b_vtgtKurama;   //!
   TBranch        *b_thetaKurama;   //!
   TBranch        *b_xtofKurama;   //!
   TBranch        *b_ytofKurama;   //!
   TBranch        *b_utofKurama;   //!
   TBranch        *b_vtofKurama;   //!
   TBranch        *b_tofsegKurama;   //!
   TBranch        *b_best_deTof;   //!
   TBranch        *b_best_TofSeg;   //!
   TBranch        *b_vpxtpc;   //!
   TBranch        *b_vpytpc;   //!
   TBranch        *b_vpztpc;   //!
   TBranch        *b_nKm;   //!
   TBranch        *b_nKp;   //!
   TBranch        *b_nKK;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vty;   //!
   TBranch        *b_vtz;   //!
   TBranch        *b_KMPX;   //!
   TBranch        *b_KMPY;   //!
   TBranch        *b_KMPZ;   //!
   TBranch        *b_KPPX;   //!
   TBranch        *b_KPPY;   //!
   TBranch        *b_KPPZ;   //!
   TBranch        *b_closeDist;   //!
   TBranch        *b_inside;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_MissMass;   //!
   TBranch        *b_MissMassCorr;   //!
   TBranch        *b_MissMassCorrDE;   //!
   TBranch        *b_MissMomx;   //!
   TBranch        *b_MissMomy;   //!
   TBranch        *b_MissMomz;   //!
   TBranch        *b_MissMomxCor;   //!
   TBranch        *b_MissMomyCor;   //!
   TBranch        *b_MissMomzCor;   //!
   TBranch        *b_thetaCM;   //!
   TBranch        *b_costCM;   //!
   TBranch        *b_Kflag;   //!
   TBranch        *b_Xiflag;   //!
   TBranch        *b_KKflag;   //!
   TBranch        *b_xkm;   //!
   TBranch        *b_ykm;   //!
   TBranch        *b_ukm;   //!
   TBranch        *b_vkm;   //!
   TBranch        *b_xkp;   //!
   TBranch        *b_ykp;   //!
   TBranch        *b_ukp;   //!
   TBranch        *b_vkp;   //!
   TBranch        *b_pOrg;   //!
   TBranch        *b_pCalc;   //!
   TBranch        *b_pCorr;   //!
   TBranch        *b_pCorrDE;   //!

   HSKKAna(TTree *tree=0);
   virtual ~HSKKAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

	void AssignTree(TTree* tree){
		TreeOut = tree;
	}
	void AssignBranches();
	void ProcessTree();
};

#endif

#ifdef HSKKAna_cxx
HSKKAna::HSKKAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DstHSKKAna.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("AllHSKKAna.root");
      }
      f->GetObject("kk",tree);

   }
   Init(tree);
}

HSKKAna::~HSKKAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HSKKAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HSKKAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HSKKAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runnum", &runnum);
   fChain->SetBranchAddress("evnum", &evnum);
   fChain->SetBranchAddress("spill", &spill);
   fChain->SetBranchAddress("trigpat", trigpat);
   fChain->SetBranchAddress("trigflag", trigflag);
   fChain->SetBranchAddress("nhBh1", &nhBh1);
   fChain->SetBranchAddress("csBh1", csBh1);
   fChain->SetBranchAddress("Bh1Seg", Bh1Seg);
   fChain->SetBranchAddress("tBh1", tBh1);
   fChain->SetBranchAddress("dtBh1", dtBh1);
   fChain->SetBranchAddress("deBh1", deBh1);
   fChain->SetBranchAddress("nhBh2", &nhBh2);
   fChain->SetBranchAddress("csBh2", csBh2);
   fChain->SetBranchAddress("Bh2Seg", Bh2Seg);
   fChain->SetBranchAddress("tBh2", tBh2);
   fChain->SetBranchAddress("t0Bh2", t0Bh2);
   fChain->SetBranchAddress("dtBh2", dtBh2);
   fChain->SetBranchAddress("deBh2", deBh2);
   fChain->SetBranchAddress("Time0Seg", &Time0Seg);
   fChain->SetBranchAddress("deTime0", &deTime0);
   fChain->SetBranchAddress("Time0", &Time0);
   fChain->SetBranchAddress("CTime0", &CTime0);
   fChain->SetBranchAddress("Btof0Seg", &Btof0Seg);
   fChain->SetBranchAddress("deBtof0", &deBtof0);
   fChain->SetBranchAddress("Btof0", &Btof0);
   fChain->SetBranchAddress("CBtof0", &CBtof0);
   fChain->SetBranchAddress("nhTof", &nhTof);
   fChain->SetBranchAddress("csTof", csTof);
   fChain->SetBranchAddress("TofSeg", TofSeg);
   fChain->SetBranchAddress("tTof", tTof);
   fChain->SetBranchAddress("dtTof", dtTof);
   fChain->SetBranchAddress("deTof", deTof);
   fChain->SetBranchAddress("nhBvh", &nhBvh);
   fChain->SetBranchAddress("BvhSeg", BvhSeg);
   fChain->SetBranchAddress("nhBft", &nhBft);
   fChain->SetBranchAddress("csBft", csBft);
   fChain->SetBranchAddress("tBft", tBft);
   fChain->SetBranchAddress("wBft", wBft);
   fChain->SetBranchAddress("BftPos", BftPos);
   fChain->SetBranchAddress("BftSeg", BftSeg);
   fChain->SetBranchAddress("nhSch", &nhSch);
   fChain->SetBranchAddress("csSch", csSch);
   fChain->SetBranchAddress("tSch", tSch);
   fChain->SetBranchAddress("wSch", wSch);
   fChain->SetBranchAddress("SchPos", SchPos);
   fChain->SetBranchAddress("SchSeg", SchSeg);
   fChain->SetBranchAddress("nlBcOut", &nlBcOut);
   fChain->SetBranchAddress("ntBcOut", &ntBcOut);
   fChain->SetBranchAddress("nhBcOut", nhBcOut);
   fChain->SetBranchAddress("chisqrBcOut", chisqrBcOut);
   fChain->SetBranchAddress("x0BcOut", x0BcOut);
   fChain->SetBranchAddress("y0BcOut", y0BcOut);
   fChain->SetBranchAddress("u0BcOut", u0BcOut);
   fChain->SetBranchAddress("v0BcOut", v0BcOut);
   fChain->SetBranchAddress("xtgtBcOut", xtgtBcOut);
   fChain->SetBranchAddress("ytgtBcOut", ytgtBcOut);
   fChain->SetBranchAddress("xbh2BcOut", xbh2BcOut);
   fChain->SetBranchAddress("ybh2BcOut", ybh2BcOut);
   fChain->SetBranchAddress("ntK18", &ntK18);
   fChain->SetBranchAddress("nhK18", nhK18);
   fChain->SetBranchAddress("chisqrK18", chisqrK18);
   fChain->SetBranchAddress("pK18", pK18);
   fChain->SetBranchAddress("xtgtHS", xtgtHS);
   fChain->SetBranchAddress("ytgtHS", ytgtHS);
   fChain->SetBranchAddress("utgtHS", utgtHS);
   fChain->SetBranchAddress("vtgtHS", vtgtHS);
   fChain->SetBranchAddress("thetaHS", thetaHS);
   fChain->SetBranchAddress("tpcHSvpx", tpcHSvpx);
   fChain->SetBranchAddress("tpcHSvpy", tpcHSvpy);
   fChain->SetBranchAddress("tpcHSvpz", tpcHSvpz);
   fChain->SetBranchAddress("nlSdcIn", &nlSdcIn);
   fChain->SetBranchAddress("ntSdcIn", &ntSdcIn);
   fChain->SetBranchAddress("nhSdcIn", nhSdcIn);
   fChain->SetBranchAddress("chisqrSdcIn", chisqrSdcIn);
   fChain->SetBranchAddress("x0SdcIn", x0SdcIn);
   fChain->SetBranchAddress("y0SdcIn", y0SdcIn);
   fChain->SetBranchAddress("u0SdcIn", u0SdcIn);
   fChain->SetBranchAddress("v0SdcIn", v0SdcIn);
   fChain->SetBranchAddress("nlSdcOut", &nlSdcOut);
   fChain->SetBranchAddress("ntSdcOut", &ntSdcOut);
   fChain->SetBranchAddress("nhSdcOut", nhSdcOut);
   fChain->SetBranchAddress("chisqrSdcOut", chisqrSdcOut);
   fChain->SetBranchAddress("x0SdcOut", x0SdcOut);
   fChain->SetBranchAddress("y0SdcOut", y0SdcOut);
   fChain->SetBranchAddress("u0SdcOut", u0SdcOut);
   fChain->SetBranchAddress("v0SdcOut", v0SdcOut);
   fChain->SetBranchAddress("ntKurama", &ntKurama);
   fChain->SetBranchAddress("nhKurama", nhKurama);
   fChain->SetBranchAddress("chisqrKurama", chisqrKurama);
   fChain->SetBranchAddress("stof", stof);
   fChain->SetBranchAddress("cstof", cstof);
   fChain->SetBranchAddress("path", path);
   fChain->SetBranchAddress("pKurama", pKurama);
   fChain->SetBranchAddress("qKurama", qKurama);
   fChain->SetBranchAddress("m2", m2);
   fChain->SetBranchAddress("xtgtKurama", xtgtKurama);
   fChain->SetBranchAddress("ytgtKurama", ytgtKurama);
   fChain->SetBranchAddress("utgtKurama", utgtKurama);
   fChain->SetBranchAddress("vtgtKurama", vtgtKurama);
   fChain->SetBranchAddress("thetaKurama", thetaKurama);
   fChain->SetBranchAddress("xtofKurama", xtofKurama);
   fChain->SetBranchAddress("ytofKurama", ytofKurama);
   fChain->SetBranchAddress("utofKurama", utofKurama);
   fChain->SetBranchAddress("vtofKurama", vtofKurama);
   fChain->SetBranchAddress("tofsegKurama", tofsegKurama);
   fChain->SetBranchAddress("best_deTof", best_deTof);
   fChain->SetBranchAddress("best_TofSeg", best_TofSeg);
   fChain->SetBranchAddress("vpxtpc", vpxtpc);
   fChain->SetBranchAddress("vpytpc", vpytpc);
   fChain->SetBranchAddress("vpztpc", vpztpc);
   fChain->SetBranchAddress("nKm", &nKm);
   fChain->SetBranchAddress("nKp", &nKp);
   fChain->SetBranchAddress("nKK", &nKK);
   fChain->SetBranchAddress("vtx", vtx);
   fChain->SetBranchAddress("vty", vty);
   fChain->SetBranchAddress("vtz", vtz);
   fChain->SetBranchAddress("KMPX", KMPX);
   fChain->SetBranchAddress("KMPY", KMPY);
   fChain->SetBranchAddress("KMPZ", KMPZ);
   fChain->SetBranchAddress("KPPX", KPPX);
   fChain->SetBranchAddress("KPPY", KPPY);
   fChain->SetBranchAddress("KPPZ", KPPZ);
   fChain->SetBranchAddress("closeDist", closeDist);
   fChain->SetBranchAddress("inside", inside);
   fChain->SetBranchAddress("theta", theta);
   fChain->SetBranchAddress("MissMass", MissMass);
   fChain->SetBranchAddress("MissMassCorr", MissMassCorr);
   fChain->SetBranchAddress("MissMassCorrDE", MissMassCorrDE);
   fChain->SetBranchAddress("MissMomx", MissMomx);
   fChain->SetBranchAddress("MissMomy", MissMomy);
   fChain->SetBranchAddress("MissMomz", MissMomz);
   fChain->SetBranchAddress("MissMomxCor", MissMomxCor);
   fChain->SetBranchAddress("MissMomyCor", MissMomyCor);
   fChain->SetBranchAddress("MissMomzCor", MissMomzCor);
   fChain->SetBranchAddress("thetaCM", thetaCM);
   fChain->SetBranchAddress("costCM", costCM);
   fChain->SetBranchAddress("Kflag", Kflag);
   fChain->SetBranchAddress("Xiflag", Xiflag);
   fChain->SetBranchAddress("KKflag", KKflag);
   fChain->SetBranchAddress("xkm", xkm);
   fChain->SetBranchAddress("ykm", ykm);
   fChain->SetBranchAddress("ukm", ukm);
   fChain->SetBranchAddress("vkm", vkm);
   fChain->SetBranchAddress("xkp", xkp);
   fChain->SetBranchAddress("ykp", ykp);
   fChain->SetBranchAddress("ukp", ukp);
   fChain->SetBranchAddress("vkp", vkp);
   fChain->SetBranchAddress("pOrg", pOrg);
   fChain->SetBranchAddress("pCalc", pCalc);
   fChain->SetBranchAddress("pCorr", pCorr);
   fChain->SetBranchAddress("pCorrDE", pCorrDE);
   Notify();
}

void HSKKAna::AssignBranches(){
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers

   TreeOut->Branch("runnum", &runnum);
   TreeOut->Branch("evnum", &evnum);
   TreeOut->Branch("spill", &spill);
   TreeOut->Branch("trigpat", trigpat);
   TreeOut->Branch("trigflag", trigflag);
   TreeOut->Branch("nhBh1", &nhBh1);
   TreeOut->Branch("csBh1", csBh1);
   TreeOut->Branch("Bh1Seg", Bh1Seg);
   TreeOut->Branch("tBh1", tBh1);
   TreeOut->Branch("dtBh1", dtBh1);
   TreeOut->Branch("deBh1", deBh1);
   TreeOut->Branch("nhBh2", &nhBh2);
   TreeOut->Branch("csBh2", csBh2);
   TreeOut->Branch("Bh2Seg", Bh2Seg);
   TreeOut->Branch("tBh2", tBh2);
   TreeOut->Branch("t0Bh2", t0Bh2);
   TreeOut->Branch("dtBh2", dtBh2);
   TreeOut->Branch("deBh2", deBh2);
   TreeOut->Branch("Time0Seg", &Time0Seg);
   TreeOut->Branch("deTime0", &deTime0);
   TreeOut->Branch("Time0", &Time0);
   TreeOut->Branch("CTime0", &CTime0);
   TreeOut->Branch("Btof0Seg", &Btof0Seg);
   TreeOut->Branch("deBtof0", &deBtof0);
   TreeOut->Branch("Btof0", &Btof0);
   TreeOut->Branch("CBtof0", &CBtof0);
   TreeOut->Branch("nhTof", &nhTof);
   TreeOut->Branch("csTof", csTof);
   TreeOut->Branch("TofSeg", TofSeg);
   TreeOut->Branch("tTof", tTof);
   TreeOut->Branch("dtTof", dtTof);
   TreeOut->Branch("deTof", deTof);
   TreeOut->Branch("nhBvh", &nhBvh);
   TreeOut->Branch("BvhSeg", BvhSeg);
   TreeOut->Branch("nhBft", &nhBft);
   TreeOut->Branch("csBft", csBft);
   TreeOut->Branch("tBft", tBft);
   TreeOut->Branch("wBft", wBft);
   TreeOut->Branch("BftPos", BftPos);
   TreeOut->Branch("BftSeg", BftSeg);
   TreeOut->Branch("nhSch", &nhSch);
   TreeOut->Branch("csSch", csSch);
   TreeOut->Branch("tSch", tSch);
   TreeOut->Branch("wSch", wSch);
   TreeOut->Branch("SchPos", SchPos);
   TreeOut->Branch("SchSeg", SchSeg);
   TreeOut->Branch("nlBcOut", &nlBcOut);
   TreeOut->Branch("ntBcOut", &ntBcOut);
   TreeOut->Branch("nhBcOut", nhBcOut);
   TreeOut->Branch("chisqrBcOut", chisqrBcOut);
   TreeOut->Branch("x0BcOut", x0BcOut);
   TreeOut->Branch("y0BcOut", y0BcOut);
   TreeOut->Branch("u0BcOut", u0BcOut);
   TreeOut->Branch("v0BcOut", v0BcOut);
   TreeOut->Branch("xtgtBcOut", xtgtBcOut);
   TreeOut->Branch("ytgtBcOut", ytgtBcOut);
   TreeOut->Branch("xbh2BcOut", xbh2BcOut);
   TreeOut->Branch("ybh2BcOut", ybh2BcOut);
   TreeOut->Branch("ntK18", &ntK18);
   TreeOut->Branch("nhK18", nhK18);
   TreeOut->Branch("chisqrK18", chisqrK18);
   TreeOut->Branch("pK18", pK18);
   TreeOut->Branch("xtgtHS", xtgtHS);
   TreeOut->Branch("ytgtHS", ytgtHS);
   TreeOut->Branch("utgtHS", utgtHS);
   TreeOut->Branch("vtgtHS", vtgtHS);
   TreeOut->Branch("thetaHS", thetaHS);
   TreeOut->Branch("tpcHSvpx", tpcHSvpx);
   TreeOut->Branch("tpcHSvpy", tpcHSvpy);
   TreeOut->Branch("tpcHSvpz", tpcHSvpz);
   TreeOut->Branch("nlSdcIn", &nlSdcIn);
   TreeOut->Branch("ntSdcIn", &ntSdcIn);
   TreeOut->Branch("nhSdcIn", nhSdcIn);
   TreeOut->Branch("chisqrSdcIn", chisqrSdcIn);
   TreeOut->Branch("x0SdcIn", x0SdcIn);
   TreeOut->Branch("y0SdcIn", y0SdcIn);
   TreeOut->Branch("u0SdcIn", u0SdcIn);
   TreeOut->Branch("v0SdcIn", v0SdcIn);
   TreeOut->Branch("nlSdcOut", &nlSdcOut);
   TreeOut->Branch("ntSdcOut", &ntSdcOut);
   TreeOut->Branch("nhSdcOut", nhSdcOut);
   TreeOut->Branch("chisqrSdcOut", chisqrSdcOut);
   TreeOut->Branch("x0SdcOut", x0SdcOut);
   TreeOut->Branch("y0SdcOut", y0SdcOut);
   TreeOut->Branch("u0SdcOut", u0SdcOut);
   TreeOut->Branch("v0SdcOut", v0SdcOut);
   TreeOut->Branch("ntKurama", &ntKurama);
   TreeOut->Branch("nhKurama", nhKurama);
   TreeOut->Branch("chisqrKurama", chisqrKurama);
   TreeOut->Branch("stof", stof);
   TreeOut->Branch("cstof", cstof);
   TreeOut->Branch("path", path);
   TreeOut->Branch("pKurama", pKurama);
   TreeOut->Branch("qKurama", qKurama);
   TreeOut->Branch("m2", m2);
   TreeOut->Branch("xtgtKurama", xtgtKurama);
   TreeOut->Branch("ytgtKurama", ytgtKurama);
   TreeOut->Branch("utgtKurama", utgtKurama);
   TreeOut->Branch("vtgtKurama", vtgtKurama);
   TreeOut->Branch("thetaKurama", thetaKurama);
   TreeOut->Branch("xtofKurama", xtofKurama);
   TreeOut->Branch("ytofKurama", ytofKurama);
   TreeOut->Branch("utofKurama", utofKurama);
   TreeOut->Branch("vtofKurama", vtofKurama);
   TreeOut->Branch("tofsegKurama", tofsegKurama);
   TreeOut->Branch("best_deTof", best_deTof);
   TreeOut->Branch("best_TofSeg", best_TofSeg);
   TreeOut->Branch("vpxtpc", vpxtpc);
   TreeOut->Branch("vpytpc", vpytpc);
   TreeOut->Branch("vpztpc", vpztpc);
   TreeOut->Branch("nKm", &nKm);
   TreeOut->Branch("nKp", &nKp);
   TreeOut->Branch("nKK", &nKK);
   TreeOut->Branch("vtx", vtx);
   TreeOut->Branch("vty", vty);
   TreeOut->Branch("vtz", vtz);
   TreeOut->Branch("KMPX", KMPX);
   TreeOut->Branch("KMPY", KMPY);
   TreeOut->Branch("KMPZ", KMPZ);
   TreeOut->Branch("KPPX", KPPX);
   TreeOut->Branch("KPPY", KPPY);
   TreeOut->Branch("KPPZ", KPPZ);
   TreeOut->Branch("closeDist", closeDist);
   TreeOut->Branch("inside", inside);
   TreeOut->Branch("theta", theta);
   TreeOut->Branch("MissMass[nKK]", MissMass);
   TreeOut->Branch("MissMassCorr", MissMassCorr);
   TreeOut->Branch("MissMassCorrDE", MissMassCorrDE);
   TreeOut->Branch("MissMomx", MissMomx);
   TreeOut->Branch("MissMomy", MissMomy);
   TreeOut->Branch("MissMomz", MissMomz);
   TreeOut->Branch("MissMomxCor", MissMomxCor);
   TreeOut->Branch("MissMomyCor", MissMomyCor);
   TreeOut->Branch("MissMomzCor", MissMomzCor);
   TreeOut->Branch("thetaCM", thetaCM);
   TreeOut->Branch("costCM", costCM);
   TreeOut->Branch("Kflag", Kflag);
   TreeOut->Branch("Xiflag[nKK]", Xiflag);
   TreeOut->Branch("KKflag[nKK]", KKflag);
   TreeOut->Branch("xkm", xkm);
   TreeOut->Branch("ykm", ykm);
   TreeOut->Branch("ukm", ukm);
   TreeOut->Branch("vkm", vkm);
   TreeOut->Branch("xkp", xkp);
   TreeOut->Branch("ykp", ykp);
   TreeOut->Branch("ukp", ukp);
   TreeOut->Branch("vkp", vkp);
   TreeOut->Branch("pOrg", pOrg);
   TreeOut->Branch("pCalc", pCalc);
   TreeOut->Branch("pCorr", pCorr);
   TreeOut->Branch("pCorrDE", pCorrDE);
}
void HSKKAna::ProcessTree(){
	for(int ient = 0;ient < fChain->GetEntries();++ient){
		fChain->GetEntry(ient);
		if(ient%10000==0){
			cout<<Form("Process : %d / %lld, ent = %lld",ient,fChain->GetEntries(),TreeOut->GetEntries())<<endl;
			
		}
		bool acpt = false;
		for(int ikk=0;ikk<nKK;++ikk){
			if(KKflag[ikk]){
				acpt = true;
				break;
			}
		}
		if(acpt)TreeOut->Fill();
	}
}

Bool_t HSKKAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HSKKAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HSKKAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HSKKAna_cxx
