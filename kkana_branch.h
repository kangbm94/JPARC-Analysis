#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

const Int_t MaxHits  = 500;
const Int_t MaxDepth = 16;
const Int_t NumOfSegBH1   = 11;
const Int_t NumOfSegBH2   =  8;
const Int_t NumOfSegBAC   =  2;
const Int_t NumOfSegSCH   = 64;
const Int_t NumOfSegTOF   = 24;
const Int_t NumOfSegHTOF  = 34;
const Int_t NumOfSegBVH   =  4;
const Int_t NumOfSegLAC   = 30;
const Int_t NumOfSegWC    = 20;
const Int_t NumOfPlaneBFT   =   2;
const Int_t NumOfSegBFT     = 160;

class kkana_branch {

 public :
  // Declaration of leaf types
  Int_t           runnum;
  Int_t           evnum;
  Int_t           spill;
  Int_t           trigpat[32];
  Int_t           trigflag[32];
  Int_t    nhBh1;
  Int_t    csBh1[NumOfSegBH1*MaxDepth];
  Double_t Bh1Seg[NumOfSegBH1*MaxDepth];
  Double_t tBh1[NumOfSegBH1*MaxDepth];
  Double_t dtBh1[NumOfSegBH1*MaxDepth];
  Double_t deBh1[NumOfSegBH1*MaxDepth];
  Int_t    nhBh2;
  Int_t    csBh2[NumOfSegBH2*MaxDepth];
  Double_t Bh2Seg[NumOfSegBH2*MaxDepth];
  Double_t tBh2[NumOfSegBH2*MaxDepth];
  Double_t t0Bh2[NumOfSegBH2*MaxDepth];
  Double_t dtBh2[NumOfSegBH2*MaxDepth];
  Double_t deBh2[NumOfSegBH2*MaxDepth];
  Double_t        Time0Seg;
  Double_t        deTime0;
  Double_t        Time0;
  Double_t        CTime0;
  Double_t        Btof0Seg;
  Double_t        deBtof0;
  Double_t        Btof0;
  Double_t        CBtof0;
  Int_t    nhTof;
  Int_t    csTof[NumOfSegTOF*MaxDepth];
  Double_t TofSeg[NumOfSegTOF*MaxDepth];
  Double_t tTof[NumOfSegTOF*MaxDepth];
  Double_t dtTof[NumOfSegTOF*MaxDepth];
  Double_t deTof[NumOfSegTOF*MaxDepth];

  Int_t           nhHtof;
  Int_t           csHtof[NumOfSegHTOF*MaxDepth];
  Double_t        HtofSeg[NumOfSegHTOF*MaxDepth];
  Double_t        tHtof[NumOfSegHTOF*MaxDepth];
  Double_t        dtHtof[NumOfSegHTOF*MaxDepth];
  Double_t        deHtof[NumOfSegHTOF*MaxDepth];
  Int_t    nhBvh;
  Double_t BvhSeg[NumOfSegBVH];

  //Fiber
  Int_t    nhBft;
  Int_t    csBft[NumOfSegBFT];
  Double_t tBft[NumOfSegBFT];
  Double_t wBft[NumOfSegBFT];
  Double_t BftPos[NumOfSegBFT];
  Double_t BftSeg[NumOfSegBFT];
  Int_t    nhSch;
  Int_t    csSch[NumOfSegSCH];
  Double_t tSch[NumOfSegSCH];
  Double_t wSch[NumOfSegSCH];
  Double_t SchPos[NumOfSegSCH];
  Double_t SchSeg[NumOfSegSCH];
  Int_t ntBcOut;
  Int_t nlBcOut;
  Int_t nhBcOut[MaxHits];
  Double_t chisqrBcOut[MaxHits];
  Double_t x0BcOut[MaxHits];
  Double_t y0BcOut[MaxHits];
  Double_t u0BcOut[MaxHits];
  Double_t v0BcOut[MaxHits];
  Double_t xtgtBcOut[MaxHits];
  Double_t ytgtBcOut[MaxHits];
  Double_t xbh2BcOut[MaxHits];
  Double_t ybh2BcOut[MaxHits];
  Int_t    ntK18;
  Int_t    nhK18[MaxHits];
  Double_t chisqrK18[MaxHits];
  Double_t pK18[MaxHits];
  Double_t xtgtK18[MaxHits];
  Double_t ytgtK18[MaxHits];
  Double_t utgtK18[MaxHits];
  Double_t vtgtK18[MaxHits];
  Double_t thetaK18[MaxHits];
  Int_t ntSdcIn;
  Int_t nlSdcIn;
  Int_t nhSdcIn[MaxHits];
  Double_t chisqrSdcIn[MaxHits];
  Double_t x0SdcIn[MaxHits];
  Double_t y0SdcIn[MaxHits];
  Double_t u0SdcIn[MaxHits];
  Double_t v0SdcIn[MaxHits];
  Int_t ntSdcOut;
  Int_t nlSdcOut;
  Int_t nhSdcOut[MaxHits];
  Double_t chisqrSdcOut[MaxHits];
  Double_t u0SdcOut[MaxHits];
  Double_t v0SdcOut[MaxHits];
  Double_t x0SdcOut[MaxHits];
  Double_t y0SdcOut[MaxHits];
  Int_t    ntKurama;
  Int_t    nhKurama[MaxHits];
  Double_t chisqrKurama[MaxHits];
  Double_t stof[MaxHits];
  Double_t cstof[MaxHits];
  Double_t path[MaxHits];
  Double_t pKurama[MaxHits];
  Double_t qKurama[MaxHits];
  Double_t m2[MaxHits];
  Double_t xtgtKurama[MaxHits];
  Double_t ytgtKurama[MaxHits];
  Double_t utgtKurama[MaxHits];
  Double_t vtgtKurama[MaxHits];
  Double_t thetaKurama[MaxHits];
  Double_t xtofKurama[MaxHits];
  Double_t ytofKurama[MaxHits];
  Double_t utofKurama[MaxHits];
  Double_t vtofKurama[MaxHits];
  Double_t tofsegKurama[MaxHits];
  Double_t best_deTof[MaxHits];
  Double_t best_TofSeg[MaxHits];
  Int_t    nKm;
  Int_t    nKp;
  Int_t    nKK;
  Double_t vtx[MaxHits];
  Double_t vty[MaxHits];
  Double_t vtz[MaxHits];
  Double_t closeDist[MaxHits];
  Int_t    inside[MaxHits];
  Double_t theta[MaxHits];
  Double_t MissMass[MaxHits];
  Double_t MissMassCorr[MaxHits];
  Double_t MissMassCorrDE[MaxHits];
  Double_t thetaCM[MaxHits];
  Double_t costCM[MaxHits];
  Int_t Kflag[MaxHits];
  Double_t xkm[MaxHits];
  Double_t ykm[MaxHits];
  Double_t ukm[MaxHits];
  Double_t vkm[MaxHits];
  Double_t xkp[MaxHits];
  Double_t ykp[MaxHits];
  Double_t ukp[MaxHits];
  Double_t vkp[MaxHits];

  Double_t pOrg[MaxHits];
  Double_t pCalc[MaxHits];
  Double_t pCorr[MaxHits];
  Double_t pCorrDE[MaxHits];

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
  TBranch        *b_nhHtof;   //!
  TBranch        *b_csHtof;   //!
  TBranch        *b_HtofSeg;   //!
  TBranch        *b_tHtof;   //!
  TBranch        *b_dtHtof;   //!
  TBranch        *b_deHtof;   //!
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
  TBranch        *b_xtgtK18;   //!
  TBranch        *b_ytgtK18;   //!
  TBranch        *b_utgtK18;   //!
  TBranch        *b_vtgtK18;   //!
  TBranch        *b_thetaK18;   //!
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
  TBranch        *b_nKm;   //!
  TBranch        *b_nKp;   //!
  TBranch        *b_nKK;   //!
  TBranch        *b_vtx;   //!
  TBranch        *b_vty;   //!
  TBranch        *b_vtz;   //!
  TBranch        *b_closeDist;   //!
  TBranch        *b_inside;   //!
  TBranch        *b_theta;   //!
  TBranch        *b_MissMass;   //!
  TBranch        *b_MissMassCorr;   //!
  TBranch        *b_MissMassCorrDE;   //!
  TBranch        *b_thetaCM;   //!
  TBranch        *b_costCM;   //!
  TBranch        *b_Kflag;   //!
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

  void Init(TTree *fChain);

};

void kkana_branch::Init(TTree *fChain)
{
  if (!fChain) return;

  fChain->SetBranchAddress("runnum", &runnum, &b_runnum);
  fChain->SetBranchAddress("evnum", &evnum, &b_evnum);
  fChain->SetBranchAddress("spill", &spill, &b_spill);
  fChain->SetBranchAddress("trigpat", trigpat, &b_trigpat);
  fChain->SetBranchAddress("trigflag", trigflag, &b_trigflag);
  fChain->SetBranchAddress("nhBh1", &nhBh1, &b_nhBh1);
  fChain->SetBranchAddress("csBh1", csBh1, &b_csBh1);
  fChain->SetBranchAddress("Bh1Seg", Bh1Seg, &b_Bh1Seg);
  fChain->SetBranchAddress("tBh1", tBh1, &b_tBh1);
  fChain->SetBranchAddress("dtBh1", dtBh1, &b_dtBh1);
  fChain->SetBranchAddress("deBh1", deBh1, &b_deBh1);
  fChain->SetBranchAddress("nhBh2", &nhBh2, &b_nhBh2);
  fChain->SetBranchAddress("csBh2", csBh2, &b_csBh2);
  fChain->SetBranchAddress("Bh2Seg", Bh2Seg, &b_Bh2Seg);
  fChain->SetBranchAddress("tBh2", tBh2, &b_tBh2);
  fChain->SetBranchAddress("t0Bh2", t0Bh2, &b_t0Bh2);
  fChain->SetBranchAddress("dtBh2", dtBh2, &b_dtBh2);
  fChain->SetBranchAddress("deBh2", deBh2, &b_deBh2);
  fChain->SetBranchAddress("Time0Seg", &Time0Seg, &b_Time0Seg);
  fChain->SetBranchAddress("deTime0", &deTime0, &b_deTime0);
  fChain->SetBranchAddress("Time0", &Time0, &b_Time0);
  fChain->SetBranchAddress("CTime0", &CTime0, &b_CTime0);
  fChain->SetBranchAddress("Btof0Seg", &Btof0Seg, &b_Btof0Seg);
  fChain->SetBranchAddress("deBtof0", &deBtof0, &b_deBtof0);
  fChain->SetBranchAddress("Btof0", &Btof0, &b_Btof0);
  fChain->SetBranchAddress("CBtof0", &CBtof0, &b_CBtof0);
  fChain->SetBranchAddress("nhTof", &nhTof, &b_nhTof);
  fChain->SetBranchAddress("csTof", csTof, &b_csTof);
  fChain->SetBranchAddress("TofSeg", TofSeg, &b_TofSeg);
  fChain->SetBranchAddress("tTof", tTof, &b_tTof);
  fChain->SetBranchAddress("dtTof", dtTof, &b_dtTof);
  fChain->SetBranchAddress("deTof", deTof, &b_deTof);
  fChain->SetBranchAddress("nhHtof", &nhHtof, &b_nhHtof);
  fChain->SetBranchAddress("csHtof", csHtof, &b_csHtof);
  fChain->SetBranchAddress("HtofSeg", HtofSeg, &b_HtofSeg);
  fChain->SetBranchAddress("tHtof", tHtof, &b_tHtof);
  fChain->SetBranchAddress("dtHtof", dtHtof, &b_dtHtof);
  fChain->SetBranchAddress("deHtof", deHtof, &b_deHtof);
  fChain->SetBranchAddress("nhBvh", &nhBvh, &b_nhBvh);
  fChain->SetBranchAddress("BvhSeg", BvhSeg, &b_BvhSeg);
  fChain->SetBranchAddress("nhBft", &nhBft, &b_nhBft);
  fChain->SetBranchAddress("csBft", csBft, &b_csBft);
  fChain->SetBranchAddress("tBft", tBft, &b_tBft);
  fChain->SetBranchAddress("wBft", wBft, &b_wBft);
  fChain->SetBranchAddress("BftPos", BftPos, &b_BftPos);
  fChain->SetBranchAddress("BftSeg", BftSeg, &b_BftSeg);
  fChain->SetBranchAddress("nhSch", &nhSch, &b_nhSch);
  fChain->SetBranchAddress("csSch", csSch, &b_csSch);
  fChain->SetBranchAddress("tSch", tSch, &b_tSch);
  fChain->SetBranchAddress("wSch", wSch, &b_wSch);
  fChain->SetBranchAddress("SchPos", SchPos, &b_SchPos);
  fChain->SetBranchAddress("SchSeg", SchSeg, &b_SchSeg);
  fChain->SetBranchAddress("nlBcOut", &nlBcOut, &b_nlBcOut);
  fChain->SetBranchAddress("ntBcOut", &ntBcOut, &b_ntBcOut);
  fChain->SetBranchAddress("nhBcOut", nhBcOut, &b_nhBcOut);
  fChain->SetBranchAddress("chisqrBcOut", chisqrBcOut, &b_chisqrBcOut);
  fChain->SetBranchAddress("x0BcOut", x0BcOut, &b_x0BcOut);
  fChain->SetBranchAddress("y0BcOut", y0BcOut, &b_y0BcOut);
  fChain->SetBranchAddress("u0BcOut", u0BcOut, &b_u0BcOut);
  fChain->SetBranchAddress("v0BcOut", v0BcOut, &b_v0BcOut);
  fChain->SetBranchAddress("xtgtBcOut", xtgtBcOut, &b_xtgtBcOut);
  fChain->SetBranchAddress("ytgtBcOut", ytgtBcOut, &b_ytgtBcOut);
  fChain->SetBranchAddress("xbh2BcOut", xbh2BcOut, &b_xbh2BcOut);
  fChain->SetBranchAddress("ybh2BcOut", ybh2BcOut, &b_ybh2BcOut);
  fChain->SetBranchAddress("ntK18", &ntK18, &b_ntK18);
  fChain->SetBranchAddress("nhK18", nhK18, &b_nhK18);
  fChain->SetBranchAddress("chisqrK18", chisqrK18, &b_chisqrK18);
  fChain->SetBranchAddress("pK18", pK18, &b_pK18);
  fChain->SetBranchAddress("xtgtK18", xtgtK18, &b_xtgtK18);
  fChain->SetBranchAddress("ytgtK18", ytgtK18, &b_ytgtK18);
  fChain->SetBranchAddress("utgtK18", utgtK18, &b_utgtK18);
  fChain->SetBranchAddress("vtgtK18", vtgtK18, &b_vtgtK18);
  fChain->SetBranchAddress("thetaK18", thetaK18, &b_thetaK18);
  fChain->SetBranchAddress("nlSdcIn", &nlSdcIn, &b_nlSdcIn);
  fChain->SetBranchAddress("ntSdcIn", &ntSdcIn, &b_ntSdcIn);
  fChain->SetBranchAddress("nhSdcIn", nhSdcIn, &b_nhSdcIn);
  fChain->SetBranchAddress("chisqrSdcIn", chisqrSdcIn, &b_chisqrSdcIn);
  fChain->SetBranchAddress("x0SdcIn", x0SdcIn, &b_x0SdcIn);
  fChain->SetBranchAddress("y0SdcIn", y0SdcIn, &b_y0SdcIn);
  fChain->SetBranchAddress("u0SdcIn", u0SdcIn, &b_u0SdcIn);
  fChain->SetBranchAddress("v0SdcIn", v0SdcIn, &b_v0SdcIn);
  fChain->SetBranchAddress("nlSdcOut", &nlSdcOut, &b_nlSdcOut);
  fChain->SetBranchAddress("ntSdcOut", &ntSdcOut, &b_ntSdcOut);
  fChain->SetBranchAddress("nhSdcOut", nhSdcOut, &b_nhSdcOut);
  fChain->SetBranchAddress("chisqrSdcOut", chisqrSdcOut, &b_chisqrSdcOut);
  fChain->SetBranchAddress("x0SdcOut", x0SdcOut, &b_x0SdcOut);
  fChain->SetBranchAddress("y0SdcOut", y0SdcOut, &b_y0SdcOut);
  fChain->SetBranchAddress("u0SdcOut", u0SdcOut, &b_u0SdcOut);
  fChain->SetBranchAddress("v0SdcOut", v0SdcOut, &b_v0SdcOut);
  fChain->SetBranchAddress("ntKurama", &ntKurama, &b_ntKurama);
  fChain->SetBranchAddress("nhKurama", nhKurama, &b_nhKurama);
  fChain->SetBranchAddress("chisqrKurama", chisqrKurama, &b_chisqrKurama);
  fChain->SetBranchAddress("stof", stof, &b_stof);
  fChain->SetBranchAddress("cstof", cstof, &b_cstof);
  fChain->SetBranchAddress("path", path, &b_path);
  fChain->SetBranchAddress("pKurama", pKurama, &b_pKurama);
  fChain->SetBranchAddress("qKurama", qKurama, &b_qKurama);
  fChain->SetBranchAddress("m2", m2, &b_m2);
  fChain->SetBranchAddress("xtgtKurama", xtgtKurama, &b_xtgtKurama);
  fChain->SetBranchAddress("ytgtKurama", ytgtKurama, &b_ytgtKurama);
  fChain->SetBranchAddress("utgtKurama", utgtKurama, &b_utgtKurama);
  fChain->SetBranchAddress("vtgtKurama", vtgtKurama, &b_vtgtKurama);
  fChain->SetBranchAddress("thetaKurama", thetaKurama, &b_thetaKurama);
  fChain->SetBranchAddress("xtofKurama", xtofKurama, &b_xtofKurama);
  fChain->SetBranchAddress("ytofKurama", ytofKurama, &b_ytofKurama);
  fChain->SetBranchAddress("utofKurama", utofKurama, &b_utofKurama);
  fChain->SetBranchAddress("vtofKurama", vtofKurama, &b_vtofKurama);
  fChain->SetBranchAddress("tofsegKurama", tofsegKurama, &b_tofsegKurama);
  fChain->SetBranchAddress("best_deTof", best_deTof, &b_best_deTof);
  fChain->SetBranchAddress("best_TofSeg", best_TofSeg, &b_best_TofSeg);
  fChain->SetBranchAddress("nKm", &nKm, &b_nKm);
  fChain->SetBranchAddress("nKp", &nKp, &b_nKp);
  fChain->SetBranchAddress("nKK", &nKK, &b_nKK);
  fChain->SetBranchAddress("vtx", vtx, &b_vtx);
  fChain->SetBranchAddress("vty", vty, &b_vty);
  fChain->SetBranchAddress("vtz", vtz, &b_vtz);
  fChain->SetBranchAddress("closeDist", closeDist, &b_closeDist);
  fChain->SetBranchAddress("inside", inside, &b_inside);
  fChain->SetBranchAddress("theta", theta, &b_theta);
  fChain->SetBranchAddress("MissMass", MissMass, &b_MissMass);
  fChain->SetBranchAddress("MissMassCorr", MissMassCorr, &b_MissMassCorr);
  fChain->SetBranchAddress("MissMassCorrDE", MissMassCorrDE, &b_MissMassCorrDE);
  fChain->SetBranchAddress("thetaCM", thetaCM, &b_thetaCM);
  fChain->SetBranchAddress("costCM", costCM, &b_costCM);
  fChain->SetBranchAddress("Kflag", Kflag, &b_Kflag);
  fChain->SetBranchAddress("xkm", xkm, &b_xkm);
  fChain->SetBranchAddress("ykm", ykm, &b_ykm);
  fChain->SetBranchAddress("ukm", ukm, &b_ukm);
  fChain->SetBranchAddress("vkm", vkm, &b_vkm);
  fChain->SetBranchAddress("xkp", xkp, &b_xkp);
  fChain->SetBranchAddress("ykp", ykp, &b_ykp);
  fChain->SetBranchAddress("ukp", ukp, &b_ukp);
  fChain->SetBranchAddress("vkp", vkp, &b_vkp);
  fChain->SetBranchAddress("pOrg", pOrg, &b_pOrg);
  fChain->SetBranchAddress("pCalc", pCalc, &b_pCalc);
  fChain->SetBranchAddress("pCorr", pCorr, &b_pCorr);
  fChain->SetBranchAddress("pCorrDE", pCorrDE, &b_pCorrDE);

}
