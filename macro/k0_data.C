#if !defined(CLING) || defined(ROOTCLING)
#include "CommonDataFormat/RangeReference.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "DataFormatsITS/TrackITS.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "DetectorsVertexing/PVertexerHelpers.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include <TLorentzVector.h>
#include "MathUtils/Cartesian.h"
#endif

using GIndex = o2::dataformats::VtxTrackIndex;
using Vec3 = ROOT::Math::SVector<double, 3>;
using V0 = o2::dataformats::V0;
using Cascade = o2::dataformats::Cascade;
using RRef = o2::dataformats::RangeReference<int, int>;
using VBracket = o2::math_utils::Bracket<int>;

using namespace o2::itsmft;

std::vector<double> selectV0Candidate(o2::vertexing::DCAFitterN<2> mFitterV0, Vec3 pV);
double getCosPA(o2::vertexing::DCAFitterN<2> mFitterV0, Vec3 pV);
bool doMCmatching(std::vector<o2::MCCompLabel> *MClabel, std::vector<std::vector<o2::MCTrack>> mcTracksMatrix, int i, int j);

void k0_data(std::string secFile = "../sim/o2trac_its.root", std::string pvFileName = "../sim/o2_primary_vertex.root")
{

  bool isMC = false;

  auto fSecondaries = TFile::Open(secFile.data());
  auto recTree = (TTree *)fSecondaries->Get("o2sim");

  std::vector<o2::MCCompLabel> *labITSvec = nullptr;
  std::vector<o2::MCTrack> *MCtracks = nullptr;

  // fill MC matrix
  std::vector<std::vector<o2::MCTrack>> mcTracksMatrix;

  std::vector<o2::its::TrackITS> *tracks = nullptr;

  recTree->SetBranchAddress("ITSTrack", &tracks);
  std::vector<o2::itsmft::ROFRecord> *rofArr = nullptr;
  recTree->SetBranchAddress("ITSTracksROF", &rofArr);

  if (isMC)
  {
    auto fMCTracks = TFile::Open("../sim/o2sim_Kine.root");
    auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
    treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
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
    recTree->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
  }

  auto pvFile = TFile::Open(pvFileName.data());
  auto pvTree = (TTree *)pvFile->Get("o2sim");

  std::vector<o2::vertexing::PVertex> *pvArr{nullptr};
  std::vector<o2::vertexing::V2TRef> *pvRefs{nullptr};
  std::vector<o2::dataformats::VtxTrackIndex> *pvIdx{nullptr};

  pvTree->SetBranchAddress("PrimaryVertex", &pvArr);
  pvTree->SetBranchAddress("PV2TrackRefs", &pvRefs);
  pvTree->SetBranchAddress("PVTrackIndices", &pvIdx);

  auto outFile = TFile("k0_inv_mass.root", "recreate");
  TH1D *histInvMass = new TH1D("k0", "k0s rec mass; M_{#pi^{+}#pi^{-}}(GeV/c^{2}); Counts", 80, 0.3, 0.8);
  TH1D *histInvMassLS = new TH1D("LS", "rec LS", 80, 0.3, 0.8);
  TH1D *cosPAH = new TH1D("cosPA", ";cos(#theta_{P})", 1000, -1, 1.);
  TH1D *cosPAHls = new TH1D("cosPAls", ";cos(#theta_{P}); Counts", 1000, -1, 1.);
  TH1D *DCAhist = new TH1D("DCA", ";DCA tracks", 1000, 0., 0.3);
  TH1D *primaryVertexZ = new TH1D("Primary Vertex Z coordinate", "; PV z coordinate (cm); Counts", 1000, 0., 70.);
  // TH1D *primaryVertexZ = new TH1D("Primary Vertex Z coordinate", "; PV z coordinate (cm); Counts", 1000, 0., 2.);

  TH1D *partPerRoFrame = new TH1D("ITS Tracks per ROframe", "ITS Tracks per ROframe; ITS Tracks per ROframe; Counts", 80, 1, 80.);
  TH1D *verticesPerRoFrame = new TH1D("Primary Vertices per ROframe", "Primary Vertices per ROframe; Primary Vertices per ROframe; Counts", 4, 1, 5.);

  cosPAHls->SetLineColor(kRed);

  o2::vertexing::DCAFitterN<2> mFitterV0;
  mFitterV0.setBz(2);

  for (int frame = 0; frame < recTree->GetEntriesFast(); frame++)
  {
    if (!recTree->GetEvent(frame) || !pvTree->GetEvent(frame))
      continue;
    if ((frame % 100) == 0)
    {
      std::cout << "Processing TF " << frame << std::endl;
    }

    std::vector<std::vector<Vec3>> pvROFs(rofArr->size());
    for (size_t idx{0}; idx < pvArr->size(); ++idx)
    {
      auto &pv{pvArr->at(idx)};
      auto &refs{pvRefs->at(idx)};
      int first{refs.getFirstEntryOfSource(0)}; // only ITS business here
      int last{first + refs.getEntriesOfSource(0) - 1};
      for (size_t rofId{0}; rofId < rofArr->size(); ++rofId)
      {
        auto &rof = rofArr->at(rofId);

        int trackStart = rof.getFirstEntry();
        int trackStop = trackStart + rof.getNEntries();
        int count{0};
        for (int i{first}; i < last; ++i)
        {
          int current{pvIdx->at(i)};
          if (current >= trackStop || current < trackStart)
          {
            count--;
          }
          else
          {
            count++;
          }
        }
        if (count > 0)
        {
          pvROFs[rofId].push_back({pv.getX(), pv.getY(), pv.getZ()});
          primaryVertexZ->Fill(pv.getZ());
        }
      }
    }

    for (size_t rofId{0}; rofId < rofArr->size(); ++rofId)
    {
      auto &rof = rofArr->at(rofId);
      int trackStart = rof.getFirstEntry();
      int trackStop = trackStart + rof.getNEntries();

      partPerRoFrame->Fill(rof.getNEntries());
      verticesPerRoFrame->Fill(pvROFs[rofId].size());

      for (int i = trackStart; i < trackStop; i++)
      {

        auto &track1 = (*tracks)[i];
        for (int j = i + 1; j < trackStop; j++)
        {
          auto &track2 = (*tracks)[j];
          if (isMC)
          {
            bool isMatched = doMCmatching(labITSvec, mcTracksMatrix, i, j);
            // if (!isMatched)
            // continue;
          }
          int nCand;
          try
          {
            nCand = mFitterV0.process(track1, track2);
          }
          catch (std::runtime_error &e)
          {
          }

          if (!nCand)
            continue;

          mFitterV0.propagateTracksToVertex();

          Vec3 selPv{0., 0., 0.};
          double prevCosPA{0.};
          double newCosPA;
          for (auto &pv : pvROFs[rofId])
          {
            newCosPA = getCosPA(mFitterV0, pv);
            if (std::abs(newCosPA) > std::abs(prevCosPA))
            {
              selPv[0] = pv[0];
              selPv[1] = pv[1];
              selPv[2] = pv[2];
            }
            prevCosPA = newCosPA;
          }

          (track1.getSign() != track2.getSign() ? cosPAH : cosPAHls)->Fill(newCosPA);

          std::vector<double> outVector = selectV0Candidate(mFitterV0, selPv);
          if (outVector[0] == -1)
          {
            continue;
          }

          double M = outVector[0];
          double cosPA = outVector[1];
          double dcaTracks = outVector[2];

          if (track1.getSign() != track2.getSign())
          {
            DCAhist->Fill(dcaTracks);
            histInvMass->Fill(M);
          }
          else
            histInvMassLS->Fill(M);
        }
      }
    }
  }

  outFile.cd();
  auto cv = TCanvas();
  cv.cd();
  histInvMass->Draw();
  histInvMassLS->SetLineColor(kRed);
  float integral = histInvMass->Integral(histInvMass->FindBin(0.4), histInvMass->FindBin(0.46));
  float integralLS = histInvMassLS->Integral(histInvMassLS->FindBin(0.4), histInvMassLS->FindBin(0.46));
  std::cout << integral / integralLS << std::endl;
  histInvMassLS->Scale(integral / integralLS);
  histInvMassLS->Draw("same");
  cv.Write();
  histInvMass->Write();
  histInvMassLS->Write();
  cosPAH->Write();
  cosPAHls->Write();
  partPerRoFrame->Write();
  verticesPerRoFrame->Write();
  DCAhist->Write();
  primaryVertexZ->Write();
  outFile.Close();
}

bool doMCmatching(std::vector<o2::MCCompLabel> *MClabel, std::vector<std::vector<o2::MCTrack>> mcTracksMatrix, int i, int j)
{
  std::vector<int> iter{i, j};
  std::vector<int> motherID;
  std::vector<int> vPDG;
  std::vector<int> vEvID;

  for (auto &k : iter)
  {

    auto lab = MClabel->at(k);

    int trackID, evID, srcID;
    bool fake;
    lab.get(trackID, evID, srcID, fake);
    if (!lab.isNoise() && lab.isValid())
    {
      motherID.push_back(mcTracksMatrix[evID][trackID].getMotherTrackId());
      vPDG.push_back(mcTracksMatrix[evID][trackID].GetPdgCode());
      vEvID.push_back(evID);
    }
    else
    {
      return false;
    }
  }
  if (motherID[0] == motherID[1] && vEvID[0] == vEvID[1] && mcTracksMatrix[vEvID[1]][motherID[1]].GetPdgCode() == 310)
  {
    return true;
  }
  return false;
}

double getCosPA(o2::vertexing::DCAFitterN<2> mFitterV0, Vec3 pV)
{

  auto sv = mFitterV0.getPCACandidate();

  auto delta = pV - sv;

  std::array<float, 3> p;
  TLorentzVector moth, prong;

  Vec3 pMom{0., 0., 0.};
  for (int i = 0; i < mFitterV0.getNProngs(); i++)
  {
    const auto &trc = mFitterV0.getTrack(i);
    trc.getPxPyPzGlo(p);
    prong.SetVectM({p[0], p[1], p[2]}, 0.139570);
    moth += prong;
    pMom[0] += p[0];
    pMom[1] += p[1];
    pMom[2] += p[2];
  }

  double cosPA{ROOT::Math::Dot(delta, pMom) / ROOT::Math::Mag(delta) / ROOT::Math::Mag(pMom)};
  return cosPA;
}

std::vector<double> selectV0Candidate(o2::vertexing::DCAFitterN<2> mFitterV0, Vec3 pV) // return invariant mass and CosPA of the selected candidate
{

  auto sv = mFitterV0.getPCACandidate();

  Vec3 delta = pV - sv;
  auto delta2 = delta[0] * delta[0] + delta[1] * delta[1];

  std::array<float, 3> p;
  TLorentzVector moth, prong;

  Vec3 pMom{0., 0., 0.};
  for (int i = 0; i < mFitterV0.getNProngs(); i++)
  {
    const auto &trc = mFitterV0.getTrack(i);
    trc.getPxPyPzGlo(p);
    prong.SetVectM({p[0], p[1], p[2]}, 0.139570);
    moth += prong;
    pMom[0] += p[0];
    pMom[1] += p[1];
    pMom[2] += p[2];
  }

  double pt2V0 = moth[0] * moth[0] + moth[1] * moth[1], prodXYv0 = delta[0] * moth[0] + delta[1] * moth[1], tDCAXY = prodXYv0 / pt2V0;
  double dcaX = delta[0] - moth[0] * tDCAXY, dcaY = delta[1] - moth[1] * tDCAXY, dca2 = dcaX * dcaX + dcaY * dcaY;

  double cosPA{ROOT::Math::Dot(delta, pMom) / ROOT::Math::Mag(delta) / ROOT::Math::Mag(pMom)};

  auto &track0 = mFitterV0.getTrack(0);
  auto &track1 = mFitterV0.getTrack(1);

  double dca_tracks = mFitterV0.getChi2AtPCACandidate();
  std::vector<double> outVec;

  if (std::abs(cosPA) < 0.99 || dca_tracks > 0.5 || delta2 < 1 || dca2 > 0.2 * 0.2 || pt2V0 < 0.1)
  {
    outVec.push_back(-1);
  }

  else if (!mFitterV0.isPropagateTracksToVertexDone() && !mFitterV0.propagateTracksToVertex())
  {
    outVec.push_back(-1);
  }

  else
  {
    outVec.push_back(moth.M());
    outVec.push_back(cosPA);
    outVec.push_back(dca_tracks);
  }

  return outVec;
}