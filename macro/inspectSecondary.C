#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "CommonDataFormat/RangeReference.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#endif

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using Cascade = o2::dataformats::Cascade;
using RRef = o2::dataformats::RangeReference<int, int>;
using VBracket = o2::math_utils::Bracket<int>;

void inspectSecondary(std::string secFile = "o2_secondary_vertex.root") {
  TH1F *histM2 = new TH1F("v0m2", "m2", 500, 0., 20.);
  auto f = TFile::Open(secFile.data());
  auto secondaries = (TTree *)f->Get("o2sim");

  std::vector<V0> *v0vec = nullptr;
  std::vector<RRef> *v02PVvec = nullptr;
  std::vector<Cascade> *cascvec = nullptr;
  std::vector<RRef> *casc2PVvec = nullptr;

  if (!secondaries->GetBranch("V0s")) {
    LOG(FATAL) << "No V0s branch!";
  }

  if (!secondaries->GetBranch("PV2V0Refs")) {
    LOG(FATAL) << "No PV2V0Refs branch!";
  }

  if (!secondaries->GetBranch("Cascades")) {
    LOG(FATAL) << "No Cascades branch!";
  }

  if (!secondaries->GetBranch("PV2CascRefs")) {
    LOG(FATAL) << "No PV2CascRefs branch!";
  }

  secondaries->SetBranchAddress("V0s", &v0vec);
  secondaries->SetBranchAddress("PV2V0Refs", &v02PVvec);
  secondaries->SetBranchAddress("Cascades", &cascvec);
  secondaries->SetBranchAddress("PV2CascRefs", &casc2PVvec);

  secondaries->GetEntry(0);

  for (auto &v0 : *v0vec) {
    // printf("%d %d %d \n", v0.posTrack().collisionId(),
    //        v0.negTrack().collisionId(), v0.collisionId());
    // v0.print();
    // histM2->Fill(v0.calcMass2(o2::track::PID::Helium3,
    // o2::track::PID::Pion));
    histM2->Fill(TMath::Sqrt(
        v0.calcMass2(o2::track::PID::Pion, o2::track::PID::Helium3)));
  }
  histM2->Draw();
}