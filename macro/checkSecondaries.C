#if !defined(CLING) || defined(ROOTCLING)
#include "CommonDataFormat/RangeReference.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
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

void checkSecondaries(std::string secFile = "o2_secondary_vertex.root") {

  // TH1F *histM2 = new TH1F("v0m2", "m2", 500, 0., 20.);
  TH1D *histPDG = new TH1D("recoPDGITS", "Reconstrcuted ITS PDG", 3, 0, 3);
  auto fMCTracks = TFile::Open("sgn_Kine.root");
  auto fSecondaries = TFile::Open(secFile.data());
  auto fITS = TFile::Open("o2trac_its.root");
  auto fITSTPC = TFile::Open("o2match_itstpc.root");
  auto fITSTPCTOF = TFile::Open("o2match_tof_itstpc.root");
  auto fTPCTOF = TFile::Open("o2match_tof_tpc.root");

  // trees
  auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
  auto treeSecondaries = (TTree *)fSecondaries->Get("o2sim");
  auto treeITS = (TTree *)fITS->Get("o2sim");

  // labels
  std::vector<o2::MCCompLabel> *labITSvec = nullptr;
  // std::vector<TrackITS> *recArr = nullptr;
  std::vector<o2::MCTrack> *MCtracks = nullptr;
  std::vector<V0> *v0vec = nullptr;
  std::vector<RRef> *v02PVvec = nullptr;
  std::vector<Cascade> *cascvec = nullptr;
  std::vector<RRef> *casc2PVvec = nullptr;

  treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
  treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
  treeSecondaries->SetBranchAddress("V0s", &v0vec);
  // treeSecondaries->SetBranchAddress("PV2V0Refs", &v02PVvec);
  // treeSecondaries->SetBranchAddress("Cascades", &cascvec);
  // treeSecondaries->SetBranchAddress("PV2CascRefs", &casc2PVvec);

  // custom containers
  std::vector<std::vector<o2::MCTrack>> mcTracksMatrix;

  auto nev = treeMCTracks->GetEntriesFast();

  mcTracksMatrix.resize(nev);
  for (int n = 0; n < nev; n++) { // loop over MC events
    treeMCTracks->GetEvent(n);
    mcTracksMatrix[n].resize(MCtracks->size());
    for (unsigned int mcI{0}; mcI < MCtracks->size(); ++mcI) {
      mcTracksMatrix[n][mcI] = MCtracks->at(mcI);
    }
  }

  for (int frame = 0; frame < treeITS->GetEntriesFast(); frame++) {
    if (!treeITS->GetEvent(frame)) {
      continue;
    }
    for (unsigned int iTrack{0}; iTrack < labITSvec->size(); ++iTrack) {
      auto lab = labITSvec->at(iTrack);
      int trackID, evID, srcID;
      bool fake;
      lab.get(trackID, evID, srcID, fake);
      if (!lab.isNoise() && lab.isValid()) {
        std::cout << mcTracksMatrix[evID][trackID].GetPdgCode() << std::endl;
        int ent = 0;
        switch (mcTracksMatrix[evID][trackID].GetPdgCode()) {
        case 1000020030:
          ent = 1;
          break;
        case 1010010030:
          ent = 2;
          break;
        }
        histPDG->Fill(ent);
      }
    }
  }
  histPDG->Draw();
  // for (auto &v0 : *v0vec)
  // {
  //     // printf("%d %d %d \n", v0.posTrack().collisionId(),
  //     //        v0.negTrack().collisionId(), v0.collisionId());
  //     // v0.print();
  //     // histM2->Fill(v0.calcMass2(o2::track::PID::Helium3,
  //     // o2::track::PID::Pion));
  //     // histM2->Fill(TMath::Sqrt(v0.calcMass2(o2::track::PID::Pion,
  //     o2::track::PID::Helium3))); std::cout <<
  //     v0.getProngID(0).getSourceName() << std::endl; std::cout <<
  //     v0.getProngID(0).getIndex() << std::endl; std::cout <<
  //     v0.getProngID(0).getSourceMask() << std::endl;
  // }
  // histM2->Draw();
}