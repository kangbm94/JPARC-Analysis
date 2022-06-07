// -*- C++ -*-

// Author: Shuhei Hayakawa

#include <fstream>
#include <iostream>
#include <sstream>

#include <TString.h>

//_____________________________________________________________________________
enum ERunInfo
  {
    kStartTime,
    kStopTime,
    kRunNumber,
    kTitle,
    kEventNumber,
    kRecorder,
    kBeamMom,
    kBeamSpill,
    kKBeamSpill,
    kPiBeamSpill,
    kBeam,
    kKBeam,
    kPiBeam,
    kSpillNumber,
    kDAQEff,
    kBeamxDAQEff,
    kDAQDuty,
    kMRPower,
    kSXDuty,
    kShsField,
    kKuramaField,
    kBLineSM1,
    kESS1_P,
    kESS1_N,
    kESS2_P,
    kESS2_N,
    // kTrigParam,
    kTpcCathodeVMon,
    kTpcGemVMon,
    kTpcGate0VMon,
    kTpcGatePVMon,
    kTpcGateMVMon,
    NRunInfo
  };

//_____________________________________________________________________________
const std::vector<TString> SRunInfo
  {
    "StartTime",
    "StopTime",
    "RunNumber",
    "RunTitle",
    "EventNumber",
    "Recorder",
    "BeamMom",
    "Beam/Spill",
    "K-Beam/Spill",
    "Pi-Beam/Spill",
    "Beam",
    "K-Beam",
    "Pi-Beam",
    "SpillNumber",
    "DAQ-Eff",
    "BeamxDAQEff",
    "DAQ-Duty",
    "MR-Power",
    "SX-Duty",
    "Shs-Field",
    "Kurama-Field",
    "BLineSM1",
    "ESS1_POS",
    "ESS1_NEG",
    "ESS2_POS",
    "ESS2_NEG",
    // "TrigParam",
    "Cathode-VMon",
    "Gem-VMon",
    "Gate0-VMon",
    "GateP-VMon",
    "GateM-VMon"
  };

//_____________________________________________________________________________
class RunSummaryMaker
{
public:
  RunSummaryMaker(const TString& file_name="RunSummary.csv",
		   const TString& delimiter=",")
    : m_file_name(file_name),
      m_file(),
      m_delimiter(delimiter)
  {
  }
  ~RunSummaryMaker(void)
  {
  }

private:
  TString                 m_file_name;
  std::ofstream           m_file;
  std::ios_base::openmode m_mode;
  TString                 m_delimiter;

public:
  Bool_t AddNewLine(std::vector<TString>& run_info);
  Bool_t MakeTitle(void);

};

//_____________________________________________________________________________
inline
Bool_t
RunSummaryMaker::AddNewLine(std::vector<TString>& run_info)
{
  m_file.open(m_file_name, std::ios::app);
  if(!m_file.is_open())
    return kFALSE;
  std::cout << "#D RunSummaryMaker::AddNewLine()" << std::endl;
  for(Int_t i=0; i<NRunInfo; ++i){
    run_info[i].ReplaceAll(",", ".");
  }
  std::stringstream ss;
  ss << run_info[kStartTime]   << m_delimiter
     << run_info[kStopTime]    << m_delimiter
     << run_info[kRunNumber]   << m_delimiter
     << run_info[kTitle]       << m_delimiter
     << run_info[kEventNumber] << m_delimiter
     << run_info[kRecorder]    << m_delimiter
     << run_info[kBeamMom]     << m_delimiter
     << run_info[kBeamSpill]   << m_delimiter
     << run_info[kKBeamSpill]  << m_delimiter
     << run_info[kPiBeamSpill] << m_delimiter
     << run_info[kBeam]        << m_delimiter
     << run_info[kKBeam]       << m_delimiter
     << run_info[kPiBeam]      << m_delimiter
     << run_info[kSpillNumber] << m_delimiter
     << run_info[kDAQEff]      << m_delimiter
     << run_info[kBeamxDAQEff] << m_delimiter
     << run_info[kDAQDuty]     << m_delimiter
     << run_info[kMRPower]     << m_delimiter
     << run_info[kSXDuty]      << m_delimiter
     << run_info[kShsField]    << m_delimiter
     << run_info[kKuramaField] << m_delimiter
     << run_info[kBLineSM1]    << m_delimiter
     << run_info[kESS1_P]      << m_delimiter
     << run_info[kESS1_N]      << m_delimiter
     << run_info[kESS2_P]      << m_delimiter
     << run_info[kESS2_N]      << m_delimiter
     // << run_info[kTrigParam]   << m_delimiter
     << run_info[kTpcCathodeVMon] << m_delimiter
     << run_info[kTpcGemVMon]     << m_delimiter
     << run_info[kTpcGate0VMon]   << m_delimiter
     << run_info[kTpcGatePVMon]   << m_delimiter
     << run_info[kTpcGateMVMon]   << m_delimiter
     << std::endl;
  std::cout << ss.str();
  m_file << ss.str();
  m_file.close();
  return kTRUE;
}

//_____________________________________________________________________________
inline
Bool_t
RunSummaryMaker::MakeTitle(void)
{
  m_file.open(m_file_name, std::ios::app);
  if(!m_file.is_open())
    return kFALSE;
  std::cout << "#D RunSummaryMaker::MakeTitle()" << std::endl;
  for(auto&& s : SRunInfo){
    m_file << s << m_delimiter;
  }
  m_file << std::endl;
  m_file.close();
  return kTRUE;
}
