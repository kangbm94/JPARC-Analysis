// -*- C++ -*-

// Author: Shuhei Hayakawa

#include "E42RunSummaryMaker.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <vector>

#include <TApplication.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TTimeStamp.h>

#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TG3DLine.h>
#include <TGNumberEntry.h>
#include <TRootEmbeddedCanvas.h>

namespace
{
const Double_t Version = 1.3;
const std::vector<TString> DataPath = {
  "/misc/raid/hddaq/e42_2021may"
};
const TString SubDataPath = "/data3/E42SubData";
const TString SavePath = SubDataPath + "/run_sheet_2021may";
const TString MtxMsTLog = SubDataPath + "/trigger_2021may/mtx.log";
const TString TriggerDir = SubDataPath + "/trigger_2021may";
const TString ScalerDir = SubDataPath + "/scaler_2021may";
const TString RunSummary = SubDataPath + "/E42RunSummary2021May.csv";
std::map<TString, Double_t> epics;
std::map<TString, TString> scaler;
std::map<TString, TString> scaler_spill;
TString run_start;
TString run_stop;
TString last_trig_param;

//_____________________________________________________________________________
// Local Functions
//_____________________________________________________________________________

//_____________________________________________________________________________
Double_t
GetEpics(const TString& key)
{
  TString cmd("caget -w 2 -t "+key);
  FILE *pipe = gSystem->OpenPipe(cmd, "r");
  if(!pipe){
    std::cerr << "#E TSystem::OpenPipe() failed" << std::endl;
    return TMath::QuietNaN();
  }
  TString line;
  line.Gets(pipe);
  gSystem->ClosePipe(pipe);
  Double_t data1, data2, data3;
  Int_t ret = std::sscanf(line.Data(), "%lf %lf %lf",
                          &data1, &data2, &data3);
  if(ret==1)
    return data1;
  else if(ret==3)
    return data2;
  else
    return TMath::QuietNaN();
}
//_____________________________________________________________________________
std::map<TString, UInt_t>
GetTriggerInfo(Int_t runnum, Bool_t wait=kTRUE)
{
  std::map<TString, UInt_t> trig;
  while(!gSystem->ProcessEvents()){
    std::ifstream ifs(TriggerDir + Form("/trigger_%05d.log", runnum));
    if(!ifs.is_open()){
      gSystem->Sleep(500);
      continue;
    }
    TString line;
    last_trig_param.Resize(0);
    while(ifs.good() && line.ReadLine(ifs)){
      if(line.IsNull() || line[0]=='#') continue;
      if(last_trig_param.IsNull())
        last_trig_param = line;
      std::istringstream iss(line.Data());
      TString buf[2];
      iss >> buf[0] >> buf[1];
      trig[buf[0]] = buf[1].Atoi();
    }
    if(trig.size() > 0)
      break;
  }
  return trig;
}

//_____________________________________________________________________________
Bool_t
IsRecorded(Int_t runnum)
{
  for(Int_t i=0, n=DataPath.size(); i<n; ++i){
    std::ifstream ifs(DataPath[i]+"/recorder.log");
    TString line;
    while(ifs.good() && line.ReadLine(ifs)){
      if(line.IsNull() || line[0] == '#')
        continue;
      std::istringstream iss(line.Data());
      TString buf[20];
      Int_t j = 0;
      while(iss >> buf[j++]);
      if(runnum == buf[1].Atoi()
         // && buf[15].Atoi() != 0
      )
        return true;
    }
  }
  return false;
}

//_____________________________________________________________________________
void
SetEpics(const TString& key)
{
  epics[key] = GetEpics(key);
}

//_____________________________________________________________________________
inline ULong64_t
Atoll(const TString& key)
{
  TString str = scaler[key];
  str.ReplaceAll(",", "");
  return str.Atoll();
}

//_____________________________________________________________________________
inline ULong64_t
DividedBySpill(const TString& key)
{
  ULong64_t val   = Atoll(key);
  ULong64_t spill = Atoll("Spill");
  if(spill!=0)
    return val/spill;
  else
    return 0.;
}

//_____________________________________________________________________________
inline TString
SeparateComma(ULong64_t num)
{
  std::vector<ULong64_t> sep_num;
  while(num/1000){
    sep_num.push_back(num%1000);
    num /= 1000;
  }

  std::stringstream ss; ss<<num;
  std::vector<ULong64_t>::reverse_iterator itr, end=sep_num.rend();
  for(itr=sep_num.rbegin(); itr!=end; ++itr){
    ss << Form(",%03llu", *itr);
  }
  return ss.str();
}

} // namespace

//_____________________________________________________________________________
class RunSheetGUI : public TGMainFrame
{
public:
  RunSheetGUI();

protected:
  // static members
  enum ETrigger     { kTRIGA, kTRIGB, kTRIGC, kTRIGD, kTRIGE, kTRIGF,
                      kCLOCK, kRESERVE2, NTrigger };
  enum EIntelligent { kMass, NIntelligent };
  enum ETarget      { kDiamond, kCH2, kEmpty, NTarget };

protected:
  TGTextButton*        fAutoButton;
  TGTextButton*        fDrawButton;
  TGTextButton*        fSaveButton;
  TGTextButton*        fPrintButton;
  TGTextButton*        fQuitButton;
  TGNumberEntry*       fRunNumber;
  TGTextEntry*         fRunTitle;
  TRootEmbeddedCanvas* fCanvas;
  TString              fSaveName;
  TLatex               fWait;
  Double_t             fWaitX;

  // Trigger
  std::vector<TString>        fTrig;
  std::vector<TGCheckButton*> fCheck;
  std::vector<TGNumberEntry*> fPS;
  // Intelligent Trigger (Matrix/Mass)
  std::vector<TString>        fIntelligent;
  std::vector<TGCheckButton*> fCheckInt;
  TGNumberEntry*              fClearPS;
  // Target
  std::vector<TString>        fTarget;
  Int_t                       fTargetID;
  TGHButtonGroup*             fTargetButtonGroup;

public:
  void   DoAuto();
  void   DoDraw();
  void   DoPrint();
  void   DoSave();
  Bool_t GetScaler(Int_t runnum);
  void   SetEnableAll();
  void   SetDisableAll(Bool_t auto_flag=kFALSE);
  void   SetRunNumber();
  void   SetRunTitle();
  void   SetTrigger();
  void   SetPreScaleFactor();
  void   SetTarget(Int_t key);

  ClassDef(RunSheetGUI, 0)
};

//_____________________________________________________________________________
RunSheetGUI::RunSheetGUI()
  : TGMainFrame(gClient->GetRoot(), 10, 10, kHorizontalFrame),
    fTrig(NTrigger),
    fCheck(NTrigger),
    fPS(NTrigger),
    fIntelligent(NIntelligent),
    fCheckInt(NIntelligent),
    fTarget(NTarget)
{
  fWait.SetNDC();
  fWait.SetTextAlign(12);
  fWait.SetTextSize(0.1);
  fWaitX = 0.280;
  SetCleanup(kDeepCleanup);
  // Controls on left
  auto controls = new TGVerticalFrame(this);
  AddFrame(controls, new TGLayoutHints(kLHintsLeft | kLHintsExpandY,
                                       5, 5, 5, 5));
  // Separator
  auto separator = new TGVertical3DLine(this);
  AddFrame(separator, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));
  // Canvas
  Double_t scale = 3.0;
  fCanvas = new TRootEmbeddedCanvas("Canvas", this, 210*scale, 297*scale);
  AddFrame(fCanvas, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY,
                                      20, 20, 20, 20));
  // Run Number
  {
    auto group = new TGGroupFrame(controls, "Run Number");
    group->SetTitlePos(TGGroupFrame::kCenter);
    fRunNumber = new TGNumberEntry(group, +0);
    group->AddFrame(fRunNumber, new TGLayoutHints(kLHintsExpandX,
                                                  0, 0, 2, 2));
    fRunNumber->Connect("ValueSet(Long_t)", "RunSheetGUI",
                        this, "SetRunNumber()");
    fRunNumber->GetNumberEntry()->Connect("ReturnPressed()", "RunSheetGUI",
                                          this, "SetRunNumber()");
    controls->AddFrame(group, new TGLayoutHints(kLHintsExpandX));
  }
  // Title
  {
    auto group = new TGGroupFrame(controls, "Run Title");
    group->SetTitlePos(TGGroupFrame::kCenter);
    fRunTitle = new TGTextEntry(group, "default title");
    fRunTitle->Resize(200);
    group->AddFrame(fRunTitle, new TGLayoutHints(kLHintsExpandX,
                                                 0, 0, 2, 2));
    fRunTitle->Connect("ReturnPressed()", "RunSheetGUI",
                       this, "SetRunTitle()");
    controls->AddFrame(group, new TGLayoutHints(kLHintsExpandX));
  }
  // Trigger
  fTrig[kTRIGA]           = "TRIG-A";
  fTrig[kTRIGB]           = "TRIG-B";
  fTrig[kTRIGC]           = "TRIG-C";
  fTrig[kTRIGD]           = "TRIG-D";
  fTrig[kTRIGE]           = "TRIG-E";
  fTrig[kTRIGF]           = "TRIG-F";
  fTrig[kCLOCK]           = "Clock";
  fTrig[kRESERVE2]        = "Reserve2";
  fIntelligent[kMass]     = "Mass";
  {
    auto group = new TGGroupFrame(controls, "Trigger");
    group->SetTitlePos(TGGroupFrame::kCenter);
    for(Int_t i=0; i<NTrigger; ++i){
      auto trigger = new TGHorizontalFrame(group);
      fCheck[i] = new TGCheckButton(trigger, fTrig[i]);
      fCheck[i]->Connect("Clicked()", "RunSheetGUI", this, "SetTrigger()");
      trigger->AddFrame(fCheck[i], new TGLayoutHints(kLHintsExpandX,
                                                     0, 0, 2, 2));
      fPS[i] = new TGNumberEntry(trigger, 1);
      fPS[i]->Connect("ValueSet(Long_t)", "RunSheetGUI",
                      this, "SetPreScaleFactor()");
      fPS[i]->GetNumberEntry()->Connect("ReturnPressed()", "RunSheetGUI",
                                        this, "SetPreScaleFactor()");
      trigger->AddFrame(fPS[i], new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));
      group->AddFrame(trigger, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));
    }
    // auto kk = new TGVButtonGroup(group, "Intelligent");
    // kk->SetTitlePos(TGGroupFrame::kCenter);
    // auto mass = new TGHorizontalFrame(kk);
    // fCheckInt[kMass] = new TGCheckButton(mass, fIntelligent[kMass]);
    // fCheckInt[kMass]->Connect("Clicked()", "RunSheetGUI",
    //                            this, "SetTrigger()");
    // fClearPS = new TGNumberEntry(mass, 1);
    // fClearPS->Connect("ValueSet(Long_t)", "RunSheetGUI",
    //                    this, "SetPreScaleFactor()");
    // fClearPS->GetNumberEntry()->Connect("ReturnPressed()", "RunSheetGUI",
    //                                      this, "SetPreScaleFactor()");
    // mass->AddFrame(fCheckInt[kMass], new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));
    // mass->AddFrame(fClearPS, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));
    // kk->AddFrame(mass, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));
    // group->AddFrame(kk, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));
    controls->AddFrame(group, new TGLayoutHints(kLHintsExpandX));
  }

  // Target
  fTarget[kDiamond] = "Diamond  ";
  fTarget[kCH2]     = "CH2  ";
  fTarget[kEmpty]   = "Empty  ";
  {
    fTargetButtonGroup = new TGHButtonGroup(controls, "Target");
    fTargetButtonGroup->SetTitlePos(TGGroupFrame::kCenter);
    for(Int_t i=0; i<NTarget; ++i){
      new TGRadioButton(fTargetButtonGroup, fTarget[i], i);
    }
    fTargetButtonGroup->Connect("Pressed(Int_t)", "RunSheetGUI", this,
                                "SetTarget(Int_t)");
    fTargetButtonGroup->SetButton(kDiamond);
    controls->AddFrame(fTargetButtonGroup, new TGLayoutHints(kLHintsExpandX));
  }

  // Button (from Bottom)
  fQuitButton = new TGTextButton(controls, "&Quit");
  controls->AddFrame(fQuitButton,
                     new TGLayoutHints(kLHintsBottom | kLHintsLeft,
                                       0, 0, 0, 5));
  fQuitButton->Connect("Pressed()", "TApplication",
                       gApplication, "Terminate()");
  fPrintButton = new TGTextButton(controls, "&Print");
  controls->AddFrame(fPrintButton,
                     new TGLayoutHints(kLHintsBottom | kLHintsExpandX,
                                       0, 0, 0, 5));
  fPrintButton->Connect("Pressed()", "RunSheetGUI", this, "DoPrint()");
  fSaveButton = new TGTextButton(controls, "&Save");
  controls->AddFrame(fSaveButton,
                     new TGLayoutHints(kLHintsBottom | kLHintsExpandX,
                                       0, 0, 0, 5));
  fSaveButton->Connect("Pressed()", "RunSheetGUI", this, "DoSave()");
  fDrawButton = new TGTextButton(controls, "&Draw");
  controls->AddFrame(fDrawButton,
                     new TGLayoutHints(kLHintsBottom | kLHintsExpandX,
                                       0, 0, 0, 5));
  fDrawButton->Connect("Pressed()", "RunSheetGUI", this, "DoDraw()");
  fAutoButton = new TGTextButton(controls, "&Auto");
  controls->AddFrame(fAutoButton,
                     new TGLayoutHints(kLHintsBottom | kLHintsExpandX,
                                       0, 0, 0, 5));
  fAutoButton->Connect("Pressed()", "RunSheetGUI", this, "DoAuto()");
  // Close Window
  Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
  DontCallClose();
  MapSubwindows();
  Resize();
  SetWMSizeHints(GetDefaultWidth(), GetDefaultHeight(), 1000, 1000, 0 ,0);
  const TString user = std::getenv("USER");
  const TString host = std::getenv("HOSTNAME");
  const TString pid  = TString::Itoa(::getpid(),10);
  const TString window =
    "Run Sheet Maker "+user+"@"+host+" (pid = "+pid+")";
  SetWindowName(window);
  MapRaised();
  // Set Default
  SetRunNumber();
  SetRunTitle();
  SetTrigger();
  SetPreScaleFactor();
}

//_____________________________________________________________________________
void
RunSheetGUI::DoAuto()
{
  if(!fAutoButton->GetString().Contains("Auto")){
    SetEnableAll();
    fAutoButton->SetText("&Auto");
    return;
  }
  std::cout << "#D RunSheetGUI::Auto()" << std::endl;
  fAutoButton->SetText("&Stop");
  SetDisableAll(true);
  Int_t runnum = fRunNumber->GetNumberEntry()->GetNumber();
  while(!gSystem->ProcessEvents()){
    SetDisableAll(true);
    run_start.Resize(0);
    run_stop.Resize(0);
    while(!gSystem->ProcessEvents()){
      if(fDrawButton->GetString().Contains("Auto"))
	return;
      for(Int_t i=0, n=DataPath.size(); i<n; ++i){
	std::ifstream ifs(DataPath[i]+"/misc/comment.txt");
	TString line;
	while(ifs.good() && line.ReadLine(ifs)){
	  if(line.IsNull() || line[0] == '#') continue;
	  std::istringstream iss(line.Data());
	  TString buf[6];
	  iss >> buf[0] >> buf[1] >> buf[2] >> buf[3] >> buf[4] >> buf[5];
	  TString num = buf[4];
	  num.ReplaceAll("]", "");
	  if(runnum != num.Atoi())
	    continue;
	  TString state = buf[5];
	  buf[1].ReplaceAll("/", "-");
	  if(state.Contains("START"))
	    run_start = buf[0] + "-" + buf[1] + " " + buf[2];
	  if(state.Contains("STOP"))
	    run_stop = buf[0] + "-" + buf[1] + " " + buf[2];
	}
      }
      if(run_start.IsNull()){
	TString p;
	for(Int_t i=0; i<5; ++i){
	  if(fAutoButton->GetString().Contains("Auto")){
	    std::cout << std::endl;
	    return;
	  }
	  p += ".";
	  std::cout << "   waiting run start " << std::setw(5) << p
		    << std::endl << "\x1b[1A";
	  gSystem->ProcessEvents();
	  gSystem->Sleep(500);
	}
      } else {
	SetRunNumber();
	break;
      }
    }
    gSystem->Sleep(10000);
    for(Int_t i=0; i<NTrigger; ++i){
      fCheck[i]->SetState(kButtonUp);
    }
    std::map<TString, UInt_t> trig = GetTriggerInfo(runnum);
    std::map<TString, UInt_t>::const_iterator itr, end;
    for(itr=trig.begin(), end=trig.end(); itr!=end; ++itr){
      TString f = itr->first;
      UInt_t  s = itr->second;
      // std::cout << f << " " << s << std::endl;
      if(f.Contains("param"))
        std::cout << "   Trigger = " << f << std::endl;
      for(Int_t i=kTRIGA; i<=kTRIGF; ++i){
        if(f == "RGN3::SEL_PSOR" && ((s>>i)&1)== 1)
          fCheck[i]->SetState(kButtonDown);
        if(f == Form("RGN3::PS_R2%c", 'A'+i))
          fPS[i]->SetNumber(s+1);
      }
    }

    if(IsRecorded(runnum)){
      DoDraw();
      gSystem->ProcessEvents();
      // DoSave();
      DoPrint();
      gSystem->ProcessEvents();
    } else {
      scaler.clear();
      scaler_spill.clear();
    }

    RunSummaryMaker run_summary(RunSummary);
    std::vector<TString> run_info(NRunInfo);
    run_info[kStartTime]   = run_start.ReplaceAll("-", "/");
    run_info[kStopTime]    = run_stop.ReplaceAll("-", "/");
    run_info[kRunNumber]   = TString::Itoa(runnum, 10);
    run_info[kTitle]       = fRunTitle->GetText();
    run_info[kEventNumber] = Form("%llu", Atoll("L2-Acc"));
    run_info[kRecorder]    = IsRecorded(runnum) ? "ON" : "OFF";
    run_info[kBeamMom]     = Form("%.2f", epics["D4:Field"]*0.299792*4.);
    run_info[kBeamSpill]   = Form("%llu", DividedBySpill("Beam"));
    run_info[kKBeamSpill]  = Form("%llu", DividedBySpill("K-Beam"));
    run_info[kPiBeamSpill] = Form("%llu", DividedBySpill("Pi-Beam"));
    run_info[kBeam]        = Form("%llu", Atoll("Beam"));
    run_info[kKBeam]       = Form("%llu", Atoll("K-Beam"));
    run_info[kPiBeam]      = Form("%llu", Atoll("Pi-Beam"));
    run_info[kSpillNumber] = Form("%llu", Atoll("Spill"));
    run_info[kDAQEff]      = scaler["DAQ-Eff"];
    run_info[kBeamxDAQEff] = Form("%.0f", Atoll("Beam") * scaler["DAQ-Eff"].Atof());
    run_info[kDAQDuty]     = scaler["Duty-Factor"];
    run_info[kMRPower]     = Form("%.2f", epics["HDSYS:MR_POWER"]);
    run_info[kSXDuty]      = Form("%.2f", epics["MRSLW:SXOPR_D2:VAL:DUTY"]);
    run_info[kShsField]    = Form("%.6f", epics["SHS:FLD:HALL"]);
    run_info[kKuramaField] = Form("%.6f", epics["KURAMA:Field"]);
    run_info[kBLineSM1]    = Form("%.0f", epics["HDPS:BSM1:SETD"]);
    run_info[kESS1_P]      = Form("%.0f", epics["HDESS:K18_ESS1:POS_VSET"]);
    run_info[kESS1_N]      = Form("%.0f", epics["HDESS:K18_ESS1:NEG_VSET"]);
    run_info[kESS2_P]      = Form("%.0f", epics["HDESS:K18_ESS2:POS_VSET"]);
    run_info[kESS2_N]      = Form("%.0f", epics["HDESS:K18_ESS2:NEG_VSET"]);
    // run_info[kTrigParam]   = gSystem->BaseName(last_trig_param);
    run_info[kTpcCathodeVMon] = Form("%.0f", epics["CAENHV:CRATE5:SLOT0:CH0:VMon"]);
    run_info[kTpcGemVMon]     = Form("%.0f", epics["CAENHV:CRATE5:SLOT0:CH1:VMon"]);
    run_info[kTpcGate0VMon]   = Form("%.0f", epics["CAENHV:CRATE5:SLOT0:CH2:VMon"]);
    run_info[kTpcGatePVMon]   = Form("%.0f", epics["CAENHV:CRATE5:SLOT0:CH3:VMon"]);
    run_info[kTpcGateMVMon]   = Form("%.0f", epics["CAENHV:CRATE5:SLOT0:CH4:VMon"]);
    // run_summary.MakeTitle();
    run_summary.AddNewLine(run_info);

    fRunNumber->SetNumber(++runnum);
    gSystem->ProcessEvents();
    SetRunNumber();
  }

  SetEnableAll();
  fAutoButton->SetText("&Auto");
}

//_____________________________________________________________________________
void
RunSheetGUI::DoDraw()
{
  if(!fDrawButton->GetString().Contains("Draw")){
    fDrawButton->SetText("&Draw");
    return;
  }

  std::cout << "#D RunSheetGUI::Draw()" << std::endl;

  fDrawButton->SetText("&Stop");

  std::vector< std::vector<TString> > mtx_mst_param;
  {
    std::ifstream ifs(MtxMsTLog);
    TString line;
    Int_t iline = 0;
    while(ifs.good() && line.ReadLine(ifs)){
      if(line.IsNull() || line[0]=='#') continue;
      std::istringstream iss(line.Data());
      std::vector<TString> buf(2);
      iss >> buf[0] >> buf[1];
      if(iline>0)
	mtx_mst_param.push_back(buf);
      iline++;
    }
  }

  SetDisableAll();
  fCanvas->GetCanvas()->Clear();
  fWait.Draw();
  fWait.SetText(fWaitX, 0.500, "Drawing.....");
  fCanvas->GetCanvas()->Modified();
  fCanvas->GetCanvas()->Update();

  const Int_t runnum = fRunNumber->GetNumberEntry()->GetNumber();

  TTimeStamp stamp;
  stamp.Add(-stamp.GetZoneOffset());

  std::cout << "   TimeStamp : " << stamp.AsString("s") << std::endl;

  ////////// Epics
  epics.clear();
  SetEpics("HDSYS:MR_POWER"); SetEpics("HDMON:MR_P3:INTENSITY");
  SetEpics("MRSLW:SXOPR_D2:VAL:DUTY"); SetEpics("MRSLW:SXOPR_D2:VAL:SpLen");
  SetEpics("D4:Field"); SetEpics("KURAMA:Field");
  SetEpics("SHS:FLD:HALL");
  SetEpics("HDPS:K18D4:SETD_MON"); SetEpics("HDPS:K18D4S:SETD_MON");
  SetEpics("HDPS:AK18D1:SETD_MON"); SetEpics("HDPS:AK11D1:SETD_MON");
  SetEpics("HDPS:BSM1:SETD");
  SetEpics("HDPS:K18CM1:SETD_MON"); SetEpics("HDPS:K18CM2:SETD_MON");
  SetEpics("HDPS:K18CM3:SETD_MON"); SetEpics("HDPS:K18CM4:SETD_MON");
  SetEpics("HDESS:K18_ESS1:NEG_VSET"); SetEpics("HDESS:K18_ESS1:POS_VSET");
  SetEpics("HDESS:K18_ESS1:NEG_IMON"); SetEpics("HDESS:K18_ESS1:POS_IMON");
  SetEpics("HDESS:K18_ESS2:NEG_VSET"); SetEpics("HDESS:K18_ESS2:POS_VSET");
  SetEpics("HDESS:K18_ESS2:NEG_IMON"); SetEpics("HDESS:K18_ESS2:POS_IMON");
  SetEpics("SLIT:IFH:LEFT"); SetEpics("SLIT:IFH:RIGHT");
  SetEpics("SLIT:IFV:UPPER"); SetEpics("SLIT:IFV:LOWER");
  SetEpics("SLIT:MOM:LEFT"); SetEpics("SLIT:MOM:RIGHT");
  SetEpics("SLIT:MASS1:UPPER"); SetEpics("SLIT:MASS1:LOWER");
  SetEpics("SLIT:MASS2:UPPER"); SetEpics("SLIT:MASS2:LOWER");
  SetEpics("AIR:BFT_TENT:TEMP"); SetEpics("AIR:BFT_TENT:HUMI");
  SetEpics("AIR:BH2_TENT:TEMP"); SetEpics("AIR:BH2_TENT:HUMI");
  SetEpics("AIR:KURAMA_TENT:TEMP"); SetEpics("AIR:KURAMA_TENT:HUMI");
  SetEpics("AIR:SCH:TEMP"); SetEpics("AIR:SCH:HUMI");
  SetEpics("CAENHV:CRATE5:SLOT0:CH0:VMon");
  SetEpics("CAENHV:CRATE5:SLOT0:CH1:VMon");
  SetEpics("CAENHV:CRATE5:SLOT0:CH2:VMon");
  SetEpics("CAENHV:CRATE5:SLOT0:CH3:VMon");
  SetEpics("CAENHV:CRATE5:SLOT0:CH4:VMon");

  if(!GetScaler(runnum)){
    std::cout << std::endl << "   stopped" << std::endl;
    fWait.SetText(0.500, 0.500, "");
    fCanvas->GetCanvas()->Modified();
    fCanvas->GetCanvas()->Update();
    fDrawButton->SetText("&Draw");
    SetEnableAll();
    return;
  }
  fCanvas->GetCanvas()->Clear();

  Double_t y = 1.0;
  TLine line;
  line.SetNDC();
  line.SetLineColor(kBlack);
  TLatex tex;
  tex.SetNDC();
  tex.SetTextAlign(12); // X:Left Y:Center
  tex.SetTextColor(kBlack);
  tex.SetTextFont(42);
  tex.SetTextSize(0.02);
  TLatex title;
  title.SetNDC();
  title.SetTextAlign(22); // X:Center Y:Center
  title.SetTextColor(kBlue+2);
  title.SetTextSize(0.02);
  //
  line.SetLineWidth(2);
  line.DrawLine(0.050, 0.050, 0.050, 0.950);
  line.DrawLine(0.050, 0.050, 0.950, 0.050);
  line.DrawLine(0.950, 0.050, 0.950, 0.950);
  line.DrawLine(0.050, 0.950, 0.950, 0.950);
  line.SetLineWidth(1);
  // title
  line.DrawLine(0.050, 0.900, 0.950, 0.900);
  line.DrawLine(0.050, 0.925, 0.200, 0.925);
  line.DrawLine(0.200, 0.050, 0.200, 0.950); //
  line.DrawLine(0.500, 0.900, 0.500, 0.950);
  title.DrawLatex(0.125, 0.9375, "E42 Run Sheet");
  title.DrawLatex(0.125, 0.9125, Form("version %.1f", Version));
  tex.SetTextColor(kRed+1);
  tex.DrawLatex(0.210, 0.9375, "Run#");
  tex.SetTextColor(kBlack);
  tex.DrawLatex(0.510, 0.9375, "Name");
  line.DrawLine(0.050, 0.860, 0.950, 0.860);
  title.DrawLatex(0.125, 0.880, "Title");
  line.DrawLine(0.050, 0.840, 0.950, 0.840);
  line.SetLineStyle(2);
  line.DrawLine(0.575, 0.775, 0.575, 0.840);
  line.DrawLine(0.660, 0.840, 0.660, 0.860);
  line.SetLineStyle(1);
  title.DrawLatex(0.125, 0.850, "Date");
  tex.DrawLatex(0.680, 0.850, "epics");
  tex.SetTextAlign(32);
  tex.DrawLatex(0.930, 0.850, stamp.AsString("s"));
  tex.SetTextAlign(12);
  tex.SetTextSize(0.05);
  tex.DrawLatex(0.300, 0.925, Form("%05d", runnum));
  tex.SetTextSize(0.025);
  tex.DrawLatex(0.220, 0.880, fRunTitle->GetText());
  tex.SetTextSize(0.020);
  // trigger
  line.DrawLine(0.050, 0.775, 0.950, 0.775);
  title.DrawLatex(0.125, 0.810, "Trigger");
  for(Int_t i=0; i<NTrigger; ++i){
    Int_t half = i/4;
    if(i == NTrigger - 1){
      if(fCheck[i]->IsOn())
        tex.DrawLatex(0.220+0.375*half+0.120, 0.830-(i-1-half*4)*0.015, "#bullet");
      else
        tex.DrawLatex(0.220+0.375*half+0.120, 0.827-(i-1-half*4)*0.015, "#circ");
    } else {
      if(fCheck[i]->IsOn())
        tex.DrawLatex(0.220+0.375*half, 0.830-(i-half*4)*0.015, "#bullet");
      else
        tex.DrawLatex(0.220+0.375*half, 0.827-(i-half*4)*0.015, "#circ");
    }
    if(i == NTrigger - 1){
      tex.DrawLatex(0.235+0.375*half+0.120, 0.830-(i-1-half*4)*0.015, fTrig[i]);
    } else {
      tex.DrawLatex(0.235+0.375*half, 0.830-(i-half*4)*0.015, fTrig[i]);
    }
    const Int_t f = fPS[i]->GetNumberEntry()->GetNumber();
    if(fCheck[i]->IsOn() && f>0)
      tex.DrawLatex(0.400+0.375*half, 0.830-(i-half*4)*0.015,
                    Form("PS  1/%d", f));
  }
  tex.DrawLatex(0.595, 0.785, gSystem->BaseName(last_trig_param));
  // tex.DrawLatex(0.595, 0.798, "#circ");
  // tex.DrawLatex(0.610, 0.800, "Others (                    )");
  // for(Int_t i=0; i<NIntelligent; ++i){
  //   if(fCheckInt[i]->IsOn())
  //     tex.DrawLatex(0.595+i*0.120, 0.785, "#bullet");
  //   else
  //     tex.DrawLatex(0.595+i*0.120, 0.783, "#circ");
  //   tex.DrawLatex(0.610+i*0.120, 0.785, fIntelligent[i]);
  // }
  // if(fCheckInt[kMass]->IsOn())
  //   tex.DrawLatex(0.790, 0.785, Form("ClearPS  1/%.0f", fClearPS->GetNumberEntry()->GetNumber()));

  // Target
  line.DrawLine(0.050, 0.755, 0.950, 0.755);
  title.DrawLatex(0.125, 0.765, "Target");
  tex.DrawLatex(0.500, 0.765, fTarget[fTargetID]);
  // Acc condition
  line.DrawLine(0.050, 0.715, 0.950, 0.715);
  title.DrawLatex(0.125, 0.735, "Acc condition");
  Double_t mr_power = epics["HDSYS:MR_POWER"];
  tex.DrawLatex(0.220, 0.742, Form("MR Power = %.2lf kW", mr_power));
  Double_t mr_inten = epics["HDMON:MR_P3:INTENSITY"];
  tex.DrawLatex(0.220, 0.728, Form("MR Intensity = %.2E PPP", mr_inten));
  Double_t sx_duty = epics["MRSLW:SXOPR_D2:VAL:DUTY"];
  tex.DrawLatex(0.490, 0.742, Form("SX Duty = %.2lf %%", sx_duty));
  Double_t spill_length = epics["MRSLW:SXOPR_D2:VAL:SpLen"];
  tex.DrawLatex(0.490, 0.728, Form("SX Spill Length = %.2lf s", spill_length));
  Double_t pk18 = epics["D4:Field"]*0.299792*4.;
  tex.DrawLatex(0.720, 0.742, Form("pK18 = %lf GeV/c", pk18));
  // Magnet
  line.DrawLine(0.050, 0.635, 0.950, 0.635);
  title.DrawLatex(0.125, 0.675, "Magnet");
  tex.DrawLatex(0.220, 0.705, "D4");
  tex.DrawLatex(0.410, 0.705, "D4s");
  tex.DrawLatex(0.600, 0.705, "D4 Field");
  tex.DrawLatex(0.220, 0.690, "KURAMA");
  tex.DrawLatex(0.600, 0.690, "KURAMA Field");
  tex.DrawLatex(0.220, 0.675, "K1.8D1");
  tex.DrawLatex(0.410, 0.675, "BLine-SM1");
  tex.DrawLatex(0.600, 0.675, "SHS Field");
  tex.DrawLatex(0.220, 0.660, "CM1");
  tex.DrawLatex(0.410, 0.660, "CM2");
  tex.DrawLatex(0.600, 0.660, "CM3");
  tex.DrawLatex(0.780, 0.660, "CM4");
  tex.DrawLatex(0.220, 0.645, "ESS1");
  tex.DrawLatex(0.600, 0.645, "ESS2");
  tex.SetTextAlign(32);
  tex.DrawLatex(0.385, 0.705, Form("%.0lf A", epics["HDPS:K18D4:SETD_MON"]));
  tex.DrawLatex(0.560, 0.705, Form("%.0lf A", epics["HDPS:K18D4S:SETD_MON"]));
  tex.DrawLatex(0.880, 0.705, Form("%lf T", epics["D4:Field"]));
  tex.DrawLatex(0.560, 0.690, Form("%.0lf A", 2400.));
  //tex.DrawLatex(0.560, 0.690, Form("%.0lf A", epics["K18:KURAMA:CSET"]));
  tex.DrawLatex(0.880, 0.690, Form("%lf T", epics["KURAMA:Field"]));
  tex.DrawLatex(0.385, 0.675, Form("%.0lf A", epics["HDPS:AK18D1:SETD_MON"]));
  tex.DrawLatex(0.560, 0.675, Form("%.0lf A", epics["HDPS:BSM1:SETD_MON"]));
  tex.DrawLatex(0.880, 0.675, Form("%lf T", epics["SHS:FLD:HALL"]));
  tex.DrawLatex(0.385, 0.660, Form("%.0lf A", epics["HDPS:K18CM1:SETD_MON"]));
  tex.DrawLatex(0.560, 0.660, Form("%.0lf A", epics["HDPS:K18CM2:SETD_MON"]));
  tex.DrawLatex(0.765, 0.660, Form("%.0lf A", epics["HDPS:K18CM3:SETD_MON"]));
  tex.DrawLatex(0.930, 0.660, Form("%.0lf A", epics["HDPS:K18CM4:SETD_MON"]));
  tex.DrawLatex(0.420, 0.645, Form("-%.0lf / %.0lf kV,",
                                   epics["HDESS:K18_ESS1:NEG_VSET"],
                                   epics["HDESS:K18_ESS1:POS_VSET"]));
  tex.DrawLatex(0.560, 0.645, Form("-%.1lf / %.1lf uA",
                                   epics["HDESS:K18_ESS1:NEG_IMON"],
                                   epics["HDESS:K18_ESS1:POS_IMON"]));
  tex.DrawLatex(0.800, 0.645, Form("-%.0lf / %.0lf kV,",
                                   epics["HDESS:K18_ESS2:NEG_VSET"],
                                   epics["HDESS:K18_ESS2:POS_VSET"]));
  tex.DrawLatex(0.930, 0.645, Form("-%.1lf / %.1lf uA",
                                   epics["HDESS:K18_ESS2:NEG_IMON"],
                                   epics["HDESS:K18_ESS2:POS_IMON"]));
  tex.SetTextAlign(12);
  // Slit
  line.DrawLine(0.050, 0.585, 0.950, 0.585);
  title.DrawLatex(0.125, 0.610, "Slit");
  tex.DrawLatex(0.220, 0.625, "IFH");
  tex.DrawLatex(0.600, 0.625, "IFV");
  tex.DrawLatex(0.220, 0.610, "MOM");
  tex.DrawLatex(0.220, 0.595, "MS1");
  tex.DrawLatex(0.600, 0.595, "MS2");
  tex.SetTextAlign(32);
  tex.DrawLatex(0.560, 0.625, Form("%.1lf mm,    %.1lf mm",
                                   epics["SLIT:IFH:LEFT"], epics["SLIT:IFH:RIGHT"]));
  tex.DrawLatex(0.900, 0.625, Form("%.1lf mm,    %.1lf mm",
                                   epics["SLIT:IFV:UPPER"], epics["SLIT:IFV:LOWER"]));
  tex.DrawLatex(0.560, 0.610, Form("%.1lf mm,    %.1lf mm",
                                   epics["SLIT:MOM:LEFT"], epics["SLIT:MOM:RIGHT"]));
  tex.DrawLatex(0.560, 0.595, Form("%.2lf mm,    %.2lf mm",
                                   epics["SLIT:MASS1:UPPER"], epics["SLIT:MASS1:LOWER"]));
  tex.DrawLatex(0.900, 0.595, Form("%.2lf mm,    %.2lf mm",
                                   epics["SLIT:MASS2:UPPER"], epics["SLIT:MASS2:LOWER"]));
  tex.SetTextAlign(12);
  // Air
  line.DrawLine(0.050, 0.560, 0.950, 0.560);
  title.DrawLatex(0.125, 0.5725, "Air");
  tex.DrawLatex(0.210, 0.5725, Form("BFT : %.1lf/%.1lf",
                                    epics["AIR:BFT_TENT:TEMP"],
                                    epics["AIR:BFT_TENT:HUMI"]));
  tex.DrawLatex(0.365, 0.5725, Form("BH2 : %.1lf/%.1lf",
                                    epics["AIR:BH2_TENT:TEMP"],
                                    epics["AIR:BH2_TENT:HUMI"]));
  tex.DrawLatex(0.520, 0.5725, Form("KURAMA : %.1lf/%.1lf",
                                    epics["AIR:KURAMA_TENT:TEMP"],
                                    epics["AIR:KURAMA_TENT:HUMI"]));
  tex.DrawLatex(0.720, 0.5725, Form("SCH : %.1lf/%.1lf",
                                    epics["AIR:SCH:TEMP"],
                                    epics["AIR:SCH:HUMI"]));
  tex.DrawLatex(0.880, 0.5725, "[#circC/%]");
  // Check
  line.DrawLine(0.050, 0.535, 0.950, 0.535);
  title.DrawLatex(0.125, 0.5475, Form("#color[%d]{Check}", kRed+1));
  tex.DrawLatex(0.220, 0.5475, "#Box Chamber HV");
  tex.DrawLatex(0.410, 0.5475, "#Box Scaler");
  tex.DrawLatex(0.530, 0.5475, "#Box Chamber efficiency");

  tex.DrawLatex(0.230, 0.850, run_start+"   -   "+run_stop);

  // Scaler /spill
  line.DrawLine(0.050, 0.475, 0.950, 0.475);
  title.DrawLatex(0.125, 0.505, "Scaler/Spill");
  tex.SetTextSize(0.022);
  tex.DrawLatex(0.220, 0.523, "BH1");
  tex.DrawLatex(0.460, 0.523, "K-Beam");
  tex.DrawLatex(0.710, 0.523, "TM");
  tex.DrawLatex(0.220, 0.505, "BH2");
  tex.DrawLatex(0.460, 0.505, "Pi-Beam");
  tex.DrawLatex(0.710, 0.505, "SY");
  tex.DrawLatex(0.220, 0.487, "L1-Req");
  tex.DrawLatex(0.460, 0.487, "L1-Acc");
  tex.DrawLatex(0.710, 0.487, "L2-Acc");
  tex.SetTextAlign(32);
  tex.DrawLatex(0.440, 0.523, scaler_spill["BH1"]);
  tex.DrawLatex(0.690, 0.523, scaler_spill["K-Beam"]);
  tex.DrawLatex(0.930, 0.523, scaler_spill["TM"]);
  tex.DrawLatex(0.440, 0.505, scaler_spill["BH2"]);
  tex.DrawLatex(0.690, 0.505, scaler_spill["Pi-Beam"]);
  tex.DrawLatex(0.930, 0.505, scaler_spill["SY"]);
  tex.DrawLatex(0.440, 0.487, scaler_spill["L1-Req"]);
  tex.DrawLatex(0.690, 0.487, scaler_spill["L1-Acc"]);
  tex.DrawLatex(0.930, 0.487, scaler_spill["L2-Acc"]);
  tex.SetTextAlign(12);
  tex.SetTextSize(0.020);
  // Scaler /run
  y = 0.295;
  line.DrawLine(0.050, y, 0.950, y);
  line.SetLineStyle(2);
  line.DrawLine(0.575, y, 0.575, y+0.180);
  line.SetLineStyle(1);
  title.DrawLatex(0.125, y+0.090, "Scaler/Run");
  tex.SetTextSize(0.022);
  {
    Int_t i=9;
    tex.DrawLatex(0.220, y+0.015+0.017*i, "Spill");
    tex.DrawLatex(0.600, y+0.015+0.017*i--, "L1-Req");
    tex.DrawLatex(0.220, y+0.015+0.017*i, "Beam");
    tex.DrawLatex(0.600, y+0.015+0.017*i--, "L1-Acc");
    tex.DrawLatex(0.220, y+0.015+0.017*i, "K-Beam");
    tex.DrawLatex(0.600, y+0.015+0.017*i--, "Clear");
    tex.DrawLatex(0.220, y+0.015+0.017*i, "Pi-Beam");
    tex.DrawLatex(0.600, y+0.015+0.017*i--, "L2-Acc");
    tex.DrawLatex(0.220, y+0.015+0.017*i, "TRIG-A");
    tex.DrawLatex(0.600, y+0.015+0.017*i--, "Beam/TM");
    tex.DrawLatex(0.220, y+0.015+0.017*i, "TRIG-B");
    tex.DrawLatex(0.600, y+0.015+0.017*i--, "L1-Req/Beam");
    tex.DrawLatex(0.220, y+0.015+0.017*i, "TRIG-C");
    tex.DrawLatex(0.600, y+0.015+0.017*i--, "Live/Real");
    tex.DrawLatex(0.220, y+0.015+0.017*i, "TRIG-D");
    tex.DrawLatex(0.600, y+0.015+0.017*i--, "DAQ-Eff");
    tex.DrawLatex(0.220, y+0.015+0.017*i, "TRIG-E");
    tex.DrawLatex(0.600, y+0.015+0.017*i--, "L2-Eff");
    tex.DrawLatex(0.220, y+0.015+0.017*i, "TRIG-F");
    tex.DrawLatex(0.600, y+0.015+0.017*i, "Duty Factor");
  }
  tex.SetTextAlign(32);
  {
    Int_t i=9;
    tex.DrawLatex(0.560, y+0.015+0.017*i, scaler["Spill"]);
    tex.DrawLatex(0.920, y+0.015+0.017*i--, scaler["L1-Req"]);
    tex.DrawLatex(0.560, y+0.015+0.017*i, scaler["Beam"]);
    tex.DrawLatex(0.920, y+0.015+0.017*i--, scaler["L1-Acc"]);
    tex.DrawLatex(0.560, y+0.015+0.017*i, scaler["K-Beam"]);
    tex.DrawLatex(0.920, y+0.015+0.017*i--, scaler["Clear"]);
    tex.DrawLatex(0.560, y+0.015+0.017*i, scaler["Pi-Beam"]);
    tex.DrawLatex(0.920, y+0.015+0.017*i--, scaler["L2-Acc"]);
    tex.DrawLatex(0.560, y+0.015+0.017*i, scaler["TRIG-A"]);
    tex.DrawLatex(0.920, y+0.015+0.017*i--, scaler["Beam/TM"]);
    tex.DrawLatex(0.560, y+0.015+0.017*i, scaler["TRIG-B"]);
    tex.DrawLatex(0.920, y+0.015+0.017*i--, scaler["L1Req/Beam"]);
    tex.DrawLatex(0.560, y+0.015+0.017*i, scaler["TRIG-C"]);
    tex.DrawLatex(0.920, y+0.015+0.017*i--, scaler["Live/Real"]);
    tex.DrawLatex(0.560, y+0.015+0.017*i, scaler["TRIG-D"]);
    tex.DrawLatex(0.920, y+0.015+0.017*i--, scaler["DAQ-Eff"]);
    tex.DrawLatex(0.560, y+0.015+0.017*i, scaler["TRIG-E"]);
    tex.DrawLatex(0.920, y+0.015+0.017*i--, scaler["L2-Eff"]);
    tex.DrawLatex(0.560, y+0.015+0.017*i, scaler["TRIG-F"]);
    tex.DrawLatex(0.920, y+0.015+0.017*i, scaler["Duty-Factor"]);
  }
  tex.SetTextAlign(12);
  tex.SetTextSize(0.020);

  // K1.8 LogBook
  y -= 0.030;
  line.DrawLine(0.050, y, 0.950, y);
  title.DrawLatex(0.125, y+0.015, "K1.8 LogBook");
  tex.DrawLatex(0.360, y+0.015, "#");
  tex.DrawLatex(0.540, y+0.015, "p");

  // Matrix/Mass trigger
  y -= 0.060;
  line.DrawLine(0.050, y, 0.950, y);
  title.DrawLatex(0.125, y+0.030, "Mtx");
  for(Int_t i=0, n=mtx_mst_param.size(); i<n; ++i){
    tex.DrawLatex(0.250, y+(3-i)*0.015, mtx_mst_param[i][0]);
    tex.DrawLatex(0.350, y+(3-i)*0.015, mtx_mst_param[i][1]);
  }

  // Comment
  y -= 0.030;
  line.DrawLine(0.050, y, 0.200, y);
  title.DrawLatex(0.125, y+0.015, "Comment");

  // // K1.8 LogBook
  // line.DrawLine(0.050, 0.125, 0.950, 0.125);
  // title.DrawLatex(0.125, 0.140, "K1.8 LogBook");
  // tex.DrawLatex(0.360, 0.140, "#");
  // tex.DrawLatex(0.540, 0.140, "p");
  // // Comment
  // line.DrawLine(0.050, 0.095, 0.200, 0.095);
  // title.DrawLatex(0.125, 0.110, "Comment");

  fCanvas->GetCanvas()->Modified();
  fCanvas->GetCanvas()->Update();
  fDrawButton->SetText("&Draw");
  SetEnableAll();
}

//_____________________________________________________________________________
void
RunSheetGUI::DoPrint()
{
  std::cout << "#D RunSheetGUI::Print()" << std::endl;
  const Int_t runnum = fRunNumber->GetNumberEntry()->GetNumber();
  fCanvas->GetCanvas()->Print(fSaveName);
  gSystem->Exec(Form("lpr -o fitplot %s", fSaveName.Data()));
}

//_____________________________________________________________________________
void
RunSheetGUI::DoSave()
{
  std::cout << "#D RunSheetGUI::Save()" << std::endl;
  const Int_t runnum = fRunNumber->GetNumberEntry()->GetNumber();
  fCanvas->GetCanvas()->Print(fSaveName);
}

//_____________________________________________________________________________
Bool_t
RunSheetGUI::GetScaler(Int_t runnum)
{
  std::cout << "#D RunSheetGUI::GetScaler()" << std::endl;

  run_start.Resize(0);
  run_stop.Resize(0);
  scaler.clear();
  scaler_spill.clear();

  if(!IsRecorded(runnum))
    return false;

  while(!gSystem->ProcessEvents()){
    if(fDrawButton->GetString().Contains("Draw"))
      return false;
    for(Int_t i=0, n=DataPath.size(); i<n; ++i){
      std::ifstream ifs(DataPath[i]+"/misc/comment.txt");
      TString line;
      while(ifs.good() && line.ReadLine(ifs)){
	if(line.IsNull() || line[0]=='#') continue;
	std::istringstream iss(line.Data());
	TString buf[6];
	iss >> buf[0] >> buf[1] >> buf[2] >> buf[3] >> buf[4] >> buf[5];

	TString num = buf[4];
	num.ReplaceAll("]","");
	if(runnum != num.Atoi())
	  continue;

	TString state = buf[5];
	buf[1].ReplaceAll("/","-");
	if(state.Contains("START"))
	  run_start = buf[0]+"-"+buf[1]+" "+buf[2];
	if(state.Contains("STOP"))
	  run_stop = buf[0]+"-"+buf[1]+" "+buf[2];
      }
    }

    if(!run_start.IsNull() && !run_stop.IsNull())
      break;

    if(run_start.IsNull()){
      std::cout << " * run is not started? : " << runnum << std::endl;
      return false;
    }

    TString p;
    for(Int_t i=0; i<5; ++i){
      if(fDrawButton->GetString().Contains("Draw"))
	return false;
      p += ".";
      std::cout << "   waiting run stop " << std::setw(5) << p
                << std::endl << "\x1b[1A";
      std::stringstream ss; ss << "Drawing" << p;
      fWait.SetText(fWaitX, 0.500, ss.str().c_str());
      fCanvas->GetCanvas()->Modified();
      fCanvas->GetCanvas()->Update();
      gSystem->ProcessEvents();
      gSystem->Sleep(500);
    }
  }

  std::cout << std::endl << " -> run stop!" << std::endl;

  TString scaler_txt;
  const Int_t max_try = 10;
  Int_t i_try = 0;
  while(!gSystem->ProcessEvents() && scaler_txt.IsNull()){
    if(fDrawButton->GetString().Contains("Draw"))
      return false;
    TString tmp = Form("%s/scaler_%05d.txt", ScalerDir.Data(), runnum);
    std::ifstream iftmp(tmp);
    if(iftmp.is_open()){
      scaler_txt = tmp;
      break;
    }
    TString p;
    for(Int_t i=0; i<5; ++i){
      if(fDrawButton->GetString().Contains("Draw"))
	return false;
      p += ".";
      std::cout << "   waiting run sync " << std::setw(5) << p
                << std::endl << "\x1b[1A";
      std::stringstream ss; ss << "Drawing" << p;
      fWait.SetText(fWaitX, 0.500, ss.str().c_str());
      fCanvas->GetCanvas()->Modified();
      fCanvas->GetCanvas()->Update();
      gSystem->ProcessEvents();
      gSystem->Sleep(500);
    }
    if(i_try == max_try) return false;
    ++i_try;
  }

  if(fDrawButton->GetString().Contains("Draw"))
    return false;
  std::cout << std::endl << " -> sync finish!" << std::endl;
  fDrawButton->SetState(kButtonDisabled);
  gSystem->ProcessEvents();
  gSystem->Sleep(2000);

  std::ifstream ifs(scaler_txt);
  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    std::cout << line << std::endl;
    std::istringstream iss(line.Data());
    TString key, val;
    iss >> key >> val;
    scaler[key] = val;
  }
  for(auto& elem: scaler){
    scaler_spill[elem.first] = SeparateComma(DividedBySpill(elem.first));
  }
  return true;
}

//_____________________________________________________________________________
void
RunSheetGUI::SetEnableAll()
{
  fAutoButton->SetState(kButtonUp);
  fDrawButton->SetState(kButtonUp);
  fSaveButton->SetState(kButtonUp);
  fPrintButton->SetState(kButtonUp);
  fRunNumber->SetState(kTRUE);
  fRunTitle->SetState(kTRUE);
  for(Int_t i=0, n=fCheck.size(); i<n; ++i){
    // fCheck[i]->SetState(kButtonUp);
    if(fCheck[i]->IsOn())
      fPS[i]->SetState(kTRUE);
  }
  // for(Int_t i=0, n=fCheckInt.size(); i<n; ++i){
  //   fCheckInt[i]->SetState(kButtonUp);
  // }
  // if(fCheckInt[kMass]->IsOn())
  //   fClearPS->SetState(kTRUE);
  fTargetButtonGroup->SetState(kTRUE);
  gSystem->ProcessEvents();
}

//_____________________________________________________________________________
void
RunSheetGUI::SetDisableAll(Bool_t auto_flag)
{
  if(auto_flag)
    fDrawButton->SetState(kButtonDisabled);
  else
    fAutoButton->SetState(kButtonDisabled);
  fSaveButton->SetState(kButtonDisabled);
  fPrintButton->SetState(kButtonDisabled);
  fRunNumber->SetState(kFALSE);
  fRunTitle->SetState(kFALSE);
  for(Int_t i=0, n=fCheck.size(); i<n; ++i){
    // fCheck[i]->SetState(kButtonDisabled);
    fPS[i]->SetState(kFALSE);
  }
  // for(Int_t i=0, n=fCheckInt.size(); i<n; ++i){
  //   fCheckInt[i]->SetState(kButtonDisabled);
  // }
  // fClearPS->SetState(kFALSE);
  fTargetButtonGroup->SetState(kFALSE);
  gSystem->ProcessEvents();
}

//_____________________________________________________________________________
void
RunSheetGUI::SetRunNumber()
{
  const Int_t runnum = fRunNumber->GetNumberEntry()->GetNumber();

  for(Int_t i=0, n=DataPath.size(); i<n; ++i){
    std::ifstream ifs(DataPath[i]+"/misc/comment.txt");
    TString line;
    while(ifs.good() && line.ReadLine(ifs)){
      if(line.IsNull() || line[0]=='#') continue;
      std::istringstream iss(line.Data());
      TString buf[5];
      iss >> buf[0] >> buf[1] >> buf[2] >> buf[3] >> buf[4];
      TString num = buf[4];
      num.ReplaceAll("]","");
      // if(num[0]=='0') num.Remove(0,1);
      line.Remove(0,40);
      if(runnum == num.Atoi()){
	fRunTitle->SetText(line);
	break;
      }
    }
  }
  std::cout << "#D RunSheetGUI::SetRunNumber()" << std::endl
	    << "   Run# " << Form("%05d", runnum) << std::endl;

  fSaveName = Form("%s/run_sheet_%05d.pdf", SavePath.Data(), runnum);
}

//_____________________________________________________________________________
void
RunSheetGUI::SetRunTitle()
{
  std::cout << "#D RunSheetGUI::SetRunTitle()" << std::endl
	    << "   " << fRunTitle->GetText() << std::endl;
}

//_____________________________________________________________________________
void
RunSheetGUI::SetTrigger()
{
  std::cout << "#D RunSheetGUI::SetTrigger()" << std::endl;
  std::cout.setf(std::ios::left);
  for(Int_t i=0; i<NTrigger; ++i){
    TString state = fCheck[i]->GetState() ? "ON" : "OFF";
    std::cout << "   " << std::setw(20) << fTrig[i]
	      << " " << state << std::endl;
    if(fCheck[i]->GetState())
      fPS[i]->SetState(kTRUE);
    else
      fPS[i]->SetState(kFALSE);
  }

  // fClearPS->SetState(fCheckInt[kMass]->GetState());

  // for(Int_t i=0; i<NIntelligent; ++i){
  //   TString state = fCheckInt[i]->GetState() ? "ON" : "OFF";
  //   std::cout << "   " << std::setw(20) << fIntelligent[i]
  //             << " " << state << std::endl;
  // }
}

//_____________________________________________________________________________
void
RunSheetGUI::SetPreScaleFactor()
{
  std::cout << "#D RunSheetGUI::SetPreScaleFactor()" << std::endl;
  std::cout.setf(std::ios::left);
  for(Int_t i=0; i<NTrigger; ++i){
    std::cout << "   " << std::setw(20) << fTrig[i]
	      << " " << fPS[i]->GetNumberEntry()->GetNumber() << std::endl;
  }
  // std::cout << "   " << std::setw(20) << "ClearPS"
  //           << " " << fClearPS->GetNumberEntry()->GetNumber() << std::endl;
}

//_____________________________________________________________________________
void
RunSheetGUI::SetTarget(Int_t key)
{
  std::cout << "#D RunSheetGUI::SetTarget()" << std::endl;
  fTargetID = key;
  std::cout << "   Target : " << fTarget[key] << std::endl;
}

//_____________________________________________________________________________
void
E42RunSheetMaker()
{
#if defined(__CINT__) || defined(__CLING__)
  std::cout << std::endl << "\033[1;33;1m"
	    << "This macro can only be executed via ACLiC"
	    << "\033[0m" << std::endl << std::endl
	    << "Usage:" << std::endl
	    << " $ root E42RunSheetMaker.C+O" << std::endl
	    << " or" << std::endl
	    << " $ root" << std::endl
	    << " [0] .x E42RunSheetMaker.C+O" << std::endl
	    << std::endl;
  gApplication->Terminate();
#endif

  new RunSheetGUI;
}
