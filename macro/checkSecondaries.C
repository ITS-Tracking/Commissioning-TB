#if !defined(CLING) || defined(ROOTCLING)
#include "CommonDataFormat/RangeReference.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#endif

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using Cascade = o2::dataformats::Cascade;
using RRef = o2::dataformats::RangeReference<int, int>;
using VBracket = o2::math_utils::Bracket<int>;

void checkSecondaries(std::string secFile = "o2_secondary_vertex.root")
{

  auto fMCTracks = TFile::Open("sgn_Kine.root");
  auto fSecondaries = TFile::Open(secFile.data());
  auto fITS = TFile::Open("o2trac_its.root");
  auto fTPC = TFile::Open("tpctracks.root");
  auto fITSTPC = TFile::Open("o2match_itstpc.root");
  auto fITSTPCTOF = TFile::Open("o2match_tof_itstpc.root");
  auto fTPCTOF = TFile::Open("o2match_tof_tpc.root");

  // trees
  auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
  auto treeSecondaries = (TTree *)fSecondaries->Get("o2sim");
  auto treeITS = (TTree *)fITS->Get("o2sim");
  auto treeTPC = (TTree *)fTPC->Get("tpcrec");
  auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");
  auto treeITSTPCTOF = (TTree *)fITSTPCTOF->Get("matchTOF");

  // labels
  std::vector<o2::MCCompLabel> *labITSvec = nullptr;
  std::vector<o2::MCCompLabel> *labTPCvec = nullptr;
  std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;
  std::vector<o2::MCCompLabel> *labITSTPCTOFvec = nullptr;

  std::vector<o2::MCTrack> *MCtracks = nullptr;
  std::vector<V0> *v0vec = nullptr;
  std::vector<RRef> *v02PVvec = nullptr;
  std::vector<Cascade> *cascvec = nullptr;
  std::vector<RRef> *casc2PVvec = nullptr;

  treeSecondaries->SetBranchAddress("V0s", &v0vec);
  treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
  treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
  treeTPC->SetBranchAddress("TPCTracksMCTruth", &labTPCvec);
  treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);
  treeITSTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTOFvec);

  std::map<std::string, std::vector<o2::MCCompLabel> *> map{{"ITS", labITSvec}, {"ITS-TPC", labITSTPCvec}, {"ITS-TPC-TOF", labITSTPCTOFvec}};

  // fill MC matrix
  std::vector<std::vector<o2::MCTrack>> mcTracksMatrix;
  auto nev = treeMCTracks->GetEntriesFast();

  mcTracksMatrix.resize(nev);
  for (int n = 0; n < nev; n++)
  { // loop over MC events
    treeMCTracks->GetEvent(n);
    mcTracksMatrix[n].resize(MCtracks->size());
    for (unsigned int mcI{0}; mcI < MCtracks->size(); ++mcI)
    {
      mcTracksMatrix[n][mcI] = MCtracks->at(mcI);
    }
  }

  auto doMatching = [&](TTree *treeDetectors,
                        std::vector<o2::MCCompLabel> *labDetectors,
                        TH1D *histo)
  {
    for (int frame = 0; frame < treeDetectors->GetEntriesFast(); frame++)
    {
      if (!treeDetectors->GetEvent(frame))
      {
        continue;
      }
      for (unsigned int iTrack{0}; iTrack < labDetectors->size(); ++iTrack)
      {
        auto lab = labDetectors->at(iTrack);
        int trackID, evID, srcID;
        bool fake;
        lab.get(trackID, evID, srcID, fake);
        if (!lab.isNoise() && lab.isValid())
        {
          int ent = 0;
          switch (mcTracksMatrix[evID][trackID].GetPdgCode())
          {
          case 1000020030:
            ent = 1;
            break;
          case 1010010030:
            ent = 2;
            break;
          }
          histo->Fill(ent);
        }
      }
    }
  };

  std::vector<TH1D *> hists(4);
  hists[0] = new TH1D("recoPDGits", "Reconstructed ITS PDG", 3, 0, 3);
  hists[1] = new TH1D("recoPDGtpc", "Reconstructed TPC PDG", 3, 0, 3);
  hists[2] = new TH1D("recoPDGitsTPC", "Reconstructed ITS+TPC PDG", 3, 0, 3);
  hists[3] = new TH1D("recoPDGitsTPCtof", "Reconstructed ITS+TPC+TOF PDG", 3, 0, 3);

  doMatching(treeITS, labITSvec, hists[0]);             // match and fill PDG histo of its tracks
  doMatching(treeTPC, labTPCvec, hists[1]);             // match and fill PDG histo of tpc tracks
  doMatching(treeITSTPC, labITSTPCvec, hists[2]);       // match and fill PDG histo of its+tpc tracks
  doMatching(treeITSTPCTOF, labITSTPCTOFvec, hists[3]); // match and fill PDG histo of its+tpc+tof tracks

  auto outFile = TFile("Secondaries.root", "recreate");
  auto *c = new TCanvas("c1", "PDG", 1000, 400);
  c->Divide(4, 1);
  for (auto iH{0}; iH < 4; ++iH)
  {
    c->cd(iH + 1);
    hists[iH]->Write();
  }
  outFile.Close();

  treeSecondaries->GetEntry();

  for (auto &v0 : *v0vec)
  {
    for (int iV0 = 0; iV0 < 2; iV0++)
    {
      if (map[v0.getProngID(iV0).getSourceName()])
      {
        auto labTrackType = map[v0.getProngID(iV0).getSourceName()];
        auto lab = labTrackType->at(v0.getProngID(iV0).getIndex());
        int trackID, evID, srcID;
        bool fake;
        lab.get(trackID, evID, srcID, fake);
        if (!lab.isNoise() && lab.isValid())
        {
          std::cout << "---------------------------------" << std::endl;
          std::cout << "Track type: " << v0.getProngID(iV0).getSourceName() << std::endl;
          std::cout << "Track PDG: " << mcTracksMatrix[evID][trackID].GetPdgCode() << std::endl;
          std::cout << "Generated MC Pt: " << mcTracksMatrix[evID][trackID].GetPt() << std::endl;
          std::cout << "Track Reco Pt: " << v0.getProng(iV0).getPt() << std::endl;
        }
      }
    }
  }
}
