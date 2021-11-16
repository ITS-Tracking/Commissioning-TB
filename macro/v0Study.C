#if !defined(CLING) || defined(ROOTCLING)
#include "CommonDataFormat/RangeReference.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/ROFRecord.h"

#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#endif

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using MCTrack = o2::MCTrack;
using Cascade = o2::dataformats::Cascade;
using RRef = o2::dataformats::RangeReference<int, int>;
using VBracket = o2::math_utils::Bracket<int>;
using namespace o2::itsmft;
using Vec3 = ROOT::Math::SVector<double, 3>;

const int motherPDG = 1010010030;
const int firstDaughterPDG = 1000020030;
const int secondDaughterPDG = -211;

// const int motherPDG = 3122;
// const int firstDaughterPDG = 2212;
// const int secondDaughterPDG = -211;

bool getITSTrack(int motherEvID, int motherTrackID, TTree *ITStree, std::vector<o2::MCCompLabel> *ITSlabel, std::vector<o2::its::TrackITS> *ITStrack);
void doMatching(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, TTree *treeDetectors, std::vector<o2::MCCompLabel> *labDetectors, TH1D *histo);
double calcMass(const V0 &v0, double dauMass[2], int dauCharges[2]);
double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG);

void v0Study()
{
    double bins[2] = {0, 0};
    double dauMass[2] = {0, 0};
    int dauCharges[2] = {2, 1};

    if (motherPDG == 1010010030)
    {
        bins[0] = 2.96;
        bins[1] = 3.04;
        dauMass[0] = 2.80839160743;
        dauMass[1] = 0.13957;
        dauCharges[0] = 2;
        dauCharges[1] = 1;
    }

    if (motherPDG == 3122)
    {
        bins[0] = 1.0;
        bins[1] = 1.2;
        dauMass[0] = 0.938272;
        dauMass[1] = 0.13957;
        dauCharges[0] = 1;
        dauCharges[1] = 1;
    }

    std::vector<TH1D *> hists(4);
    hists[0] = new TH1D("recoPDGits", "Reconstructed ITS PDG;;Efficiency", 3, 0, 3);
    hists[1] = new TH1D("recoPDGitsTPC", "Reconstructed ITS-TPC PDG;;Efficiency", 3, 0, 3);
    hists[2] = new TH1D("recoPDGtpcTOF", "Reconstructed TPC-TOF PDG;;Efficiency", 3, 0, 3);
    hists[3] = new TH1D("recoPDGitsTPCTOF", "Reconstructed ITS-TPC-TOF PDG;;Efficiency", 3, 0, 3);
    TH1D *histInvMass = new TH1D("V0 invariant mass", "; V0 Mass (GeV/c^{2}); Counts", 30, bins[0], bins[1]);
    TH1D *histV0radius = new TH1D("V0 radius", "; V0 Radius (cm); Counts", 30, 0, 40);

    auto fMCTracks = TFile::Open("o2sim_Kine.root");
    auto fSecondaries = TFile::Open("o2_secondary_vertex.root");
    auto fITS = TFile::Open("o2trac_its.root");
    auto fITSTPC = TFile::Open("o2match_itstpc.root");
    auto fTPCTOF = TFile::Open("o2match_tof_tpc.root");
    auto fITSTPCTOF = TFile::Open("o2match_tof_itstpc.root");

    // Trees
    auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
    auto treeSecondaries = (TTree *)fSecondaries->Get("o2sim");
    auto treeITS = (TTree *)fITS->Get("o2sim");
    auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");
    auto treeITSTPCTOF = (TTree *)fITSTPCTOF->Get("matchTOF");
    auto treeTPCTOF = (TTree *)fTPCTOF->Get("matchTOF");

    // Tracks
    std::vector<o2::MCTrack> *MCtracks = nullptr;
    std::vector<V0> *v0vec = nullptr;
    std::vector<o2::its::TrackITS> *ITStracks = nullptr;
    std::vector<o2::itsmft::ROFRecord> *rofArr = nullptr;

    // Labels
    std::vector<o2::MCCompLabel> *labITSvec = nullptr;
    std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;
    std::vector<o2::MCCompLabel> *labITSTPCTOFvec = nullptr;
    std::vector<o2::MCCompLabel> *labTPCTOFvec = nullptr;

    treeSecondaries->SetBranchAddress("V0s", &v0vec);
    treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
    treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);
    treeTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTOFvec);
    treeITSTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTOFvec);
    treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
    treeITS->SetBranchAddress("ITSTrack", &ITStracks);
    treeITS->SetBranchAddress("ITSTracksROF", &rofArr);

    std::map<std::string, std::vector<o2::MCCompLabel> *> map{{"ITS", labITSvec}, {"ITS-TPC", labITSTPCvec}, {"TPC-TOF", labTPCTOFvec}, {"ITS-TPC-TOF", labITSTPCTOFvec}};

    // fill MC matrix
    int injectedParticles = 0;
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
            if (MCtracks->at(mcI).GetPdgCode() == motherPDG)
            {
                injectedParticles++;
            }
        }
    }

    doMatching(mcTracksMatrix, treeITS, labITSvec, hists[0]);
    doMatching(mcTracksMatrix, treeITSTPC, labITSTPCvec, hists[1]);
    doMatching(mcTracksMatrix, treeTPCTOF, labTPCTOFvec, hists[2]);
    doMatching(mcTracksMatrix, treeITSTPCTOF, labITSTPCTOFvec, hists[3]);

    auto outFile = TFile("v0study.root", "recreate");

    treeSecondaries->GetEntry();
    treeITS->GetEntry();
    treeITSTPC->GetEntry();
    treeTPCTOF->GetEntry();
    treeITSTPCTOF->GetEntry();

    int counter = 0;
    for (auto &v0 : *v0vec)
    {
        std::vector<int> motherIDvec;
        std::vector<int> daughterIDvec;
        std::vector<int> evIDvec;

        for (int iV0 = 0; iV0 < 2; iV0++)
        {
            if (map[v0.getProngID(iV0).getSourceName()])
            {

                auto labTrackType = map[v0.getProngID(iV0).getSourceName()];
                auto lab = labTrackType->at(v0.getProngID(iV0).getIndex());

                int trackID, evID, srcID;
                bool fake;
                lab.get(trackID, evID, srcID, fake);
                if (!lab.isNoise() && lab.isValid() && lab.isCorrect())
                {
                    auto motherID = mcTracksMatrix[evID][trackID].getMotherTrackId();
                    motherIDvec.push_back(motherID);
                    daughterIDvec.push_back(trackID);
                    evIDvec.push_back(evID);
                }
            }
        }

        if (motherIDvec.size() < 2)
            continue;
        if (motherIDvec[0] != motherIDvec[1] || evIDvec[0] != evIDvec[1])
            continue;
        if (motherIDvec[0] <= 0 || motherIDvec[0] > 10000)
            continue;

        int pdg0 = mcTracksMatrix[evIDvec[0]][daughterIDvec[0]].GetPdgCode();
        int pdg1 = mcTracksMatrix[evIDvec[0]][daughterIDvec[1]].GetPdgCode();

        histInvMass->Fill(calcMass(v0, dauMass, dauCharges));
        histV0radius->Fill(TMath::Sqrt(v0.calcR2()));

        if (pdg0 != firstDaughterPDG && pdg0 != secondDaughterPDG)
            continue;
        if (pdg1 != firstDaughterPDG && pdg1 != secondDaughterPDG)
            continue;
        std::cout << "---------------------------------" << std::endl;

        counter++;
        std::cout << evIDvec[0] << ", " << motherIDvec[0] << ", " << motherIDvec[1] << std::endl;
        std::cout << "Common mother found, PDG: " << mcTracksMatrix[evIDvec[0]][motherIDvec[0]].GetPdgCode() << std::endl;
        std::cout << "Daughter 0, PDG: " << pdg0 << ", Pt: " << mcTracksMatrix[evIDvec[0]][daughterIDvec[0]].GetPt() << std::endl;
        std::cout << "Daughter 0, Rec Pt: " << v0.getProng(0).getPt() << ", Track type: " << v0.getProngID(0).getSourceName() << std::endl;

        std::cout << "Daughter 1, PDG: " << pdg1 << ", Pt: " << mcTracksMatrix[evIDvec[0]][daughterIDvec[1]].GetPt() << std::endl;
        std::cout << "Daughter 1, Rec Pt: " << v0.getProng(1).getPt() << ", Track type: " << v0.getProngID(1).getSourceName() << std::endl;

        getITSTrack(evIDvec[0], motherIDvec[0], treeITS, labITSvec, ITStracks);
        auto motherTrack = mcTracksMatrix[evIDvec[0]][motherIDvec[0]];
        std::cout << "Did ITS see mother track (hits)? : " << motherTrack.leftTrace(0) << std::endl;
        std::cout << "Counter: " << counter << std::endl;
    }

    const char *labels[3] = {std::to_string(motherPDG).data(), std::to_string(firstDaughterPDG).data(), std::to_string(secondDaughterPDG).data()};
    double norm = 1 / double(injectedParticles);
    for (auto iH{0}; iH < 4; ++iH)
    {
        for (auto iLab{0}; iLab < 3; ++iLab)
        {
            hists[iH]->GetXaxis()->SetBinLabel(iLab + 1, labels[iLab]);
        }
        hists[iH]->Scale(norm);
        hists[iH]->Write();
    }
    histInvMass->Write();
    histV0radius->Write();
    outFile.Close();
}

bool getITSTrack(int motherEvID, int motherTrackID, TTree *ITStree, std::vector<o2::MCCompLabel> *ITSlabel, std::vector<o2::its::TrackITS> *ITStrack)
{
    // o2::its::TrackITS *motherTrack{nullptr};
    for (int frame = 0; frame < ITStree->GetEntriesFast(); frame++)
    {
        if (!ITStree->GetEvent(frame) || !ITStree->GetEvent(frame))
            continue;
        if (!ITStree->GetEvent(frame))
        {
            continue;
        }
        for (unsigned int iTrack{0}; iTrack < ITSlabel->size(); ++iTrack)
        {
            auto lab = ITSlabel->at(iTrack);
            int trackID, evID, srcID;
            bool fake;
            lab.get(trackID, evID, srcID, fake);
            if (!lab.isNoise() && lab.isValid())
            {
                if (evID == motherEvID and trackID == motherTrackID)
                {
                    auto motherTrack = ITStrack->at(iTrack);
                    std::cout << "ITS track found! " << std::endl;
                }
            }
        }
    }
    return 0;
};

void doMatching(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, TTree *treeDetectors,
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
                int ent = -1;
                switch (mcTracksMatrix[evID][trackID].GetPdgCode())
                {
                case secondDaughterPDG:
                    ent = 2;
                    break;
                case firstDaughterPDG:
                    ent = 1;
                    break;
                case motherPDG:
                    ent = 0;
                    break;
                }
                histo->Fill(ent);
            }
        }
    }
};

double calcMass(const V0 &v0, double dauMass[2], int dauCharges[2])
{
    std::vector<o2::dataformats::V0::Track> dauTracks = {v0.getProng(0), v0.getProng(1)};
    TLorentzVector moth, prong;
    std::array<float, 3> p;
    for (int i = 0; i < 2; i++)
    {
        auto &track = dauTracks[i];
        auto &mass = dauMass[i];
        track.getPxPyPzGlo(p);
        int charge = dauCharges[i];
        prong.SetVectM({charge * p[0], charge * p[1], charge * p[2]}, mass);
        moth += prong;
    }
    return moth.M();
}

double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
{
    auto idStart = motherTrack.getFirstDaughterTrackId();
    auto idStop = motherTrack.getLastDaughterTrackId();
    for (auto iD{idStart}; iD < idStop; ++iD)
    {
        auto dauTrack = MCTracks->at(iD);
        if (dauTrack.GetPdgCode() == dauPDG)
        {
            auto decLength = TMath::Sqrt((dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) * (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) + (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) * (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) + (dauTrack.GetStartVertexCoordinatesZ() - motherTrack.GetStartVertexCoordinatesZ()) * (dauTrack.GetStartVertexCoordinatesZ() - motherTrack.GetStartVertexCoordinatesZ()));
            return decLength;
        }
    }
    return -1;
}