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
#include "TLegend.h"
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

double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG);

void decLength()
{

    TH1D *histDecLength = new TH1D("Generated V0 decay length", "Generated Dec Length; Decay radius (cm); Efficiency", 50, 0, 80);
    TH1D *histDecLengthRec = new TH1D("ITS Reco V0 decay length", "ITS Reco Dec Length; Decay radius (cm); Efficiency", 50, 0, 80);
    TH1D *histDecLengthRecCorrect = new TH1D("ITS Reco V0 decay length no fake", "ITS Reco Dec Length no fakee; Decay radius; Efficiency", 50, 0, 80);
    TH1D *histDecLengthHits = new TH1D("ITS Hits V0 decay length", "ITS Hits Dec Length; ITS Hits Decay Length (cm); Efficiency", 50, 0, 80);
    TH1D *histFakeHits = new TH1D("ITS Fake hits", "ITS fake hits; Layer Number; Fake hits", 7, -0.5, 6.5);

    auto fMCTracks = TFile::Open("o2sim_Kine.root");
    auto fSecondaries = TFile::Open("o2_secondary_vertex.root");
    auto fITS = TFile::Open("o2trac_its.root");

    // Trees
    auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
    auto treeSecondaries = (TTree *)fSecondaries->Get("o2sim");
    auto treeITS = (TTree *)fITS->Get("o2sim");

    // Tracks
    std::vector<o2::MCTrack> *MCtracks = nullptr;
    std::vector<V0> *v0vec = nullptr;
    std::vector<o2::its::TrackITS> *ITStracks = nullptr;
    std::vector<o2::itsmft::ROFRecord> *rofArr = nullptr;

    // Labels
    std::vector<o2::MCCompLabel> *labITSvec = nullptr;

    treeSecondaries->SetBranchAddress("V0s", &v0vec);
    treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
    treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
    treeITS->SetBranchAddress("ITSTrack", &ITStracks);
    treeITS->SetBranchAddress("ITSTracksROF", &rofArr);

    std::map<std::string, std::vector<o2::MCCompLabel> *> map{{"ITS", labITSvec}};

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

                histDecLength->Fill(calcDecLength(MCtracks, MCtracks->at(mcI), firstDaughterPDG));
                if (mcTracksMatrix[n][mcI].leftTrace(0))
                    histDecLengthHits->Fill(calcDecLength(MCtracks, MCtracks->at(mcI), firstDaughterPDG));
            }
        }
    }

    // o2::its::TrackITS *motherTrack{nullptr};
    for (int frame = 0; frame < treeITS->GetEntriesFast(); frame++)
    {
        if (!treeITS->GetEvent(frame))
        {
            continue;
        }
        for (unsigned int iTrack{0}; iTrack < labITSvec->size(); ++iTrack)
        {
            auto lab = labITSvec->at(iTrack);
            int trackID, evID, srcID;
            bool fake;
            lab.get(trackID, evID, srcID, fake);
            if (!lab.isNoise() && lab.isValid())
            {
                auto mcTrack = mcTracksMatrix[evID][trackID];
                if (mcTrack.GetPdgCode() == motherPDG)
                {
                    histDecLengthRec->Fill(calcDecLength(&mcTracksMatrix[evID], mcTrack, firstDaughterPDG));
                    if (lab.isCorrect())
                        histDecLengthRecCorrect->Fill(calcDecLength(&mcTracksMatrix[evID], mcTrack, firstDaughterPDG));

                    auto hypITSTrack = ITStracks->at(iTrack);
                    for (auto i{0}; i < 7; i++)
                    {
                        if (hypITSTrack.isFakeOnLayer(i))
                            histFakeHits->Fill(i);
                    }
                }
            }
        }
    }

    histDecLengthRec->Divide(histDecLength);
    histDecLengthRecCorrect->Divide(histDecLength);
    histDecLengthHits->Divide(histDecLength);



    auto outFile = TFile("decLength.root", "recreate");
    auto *c = new TCanvas("c1", "decLength", 1000, 400);
    histDecLengthRec->SetStats(0);
    histDecLengthRec->SetLineColor(kRed);
    histDecLengthRec->SetLineWidth(2);

    histDecLengthRecCorrect->SetLineColor(kOrange);
    histDecLengthRecCorrect->SetLineWidth(2);

    histDecLengthHits->SetLineColor(kGreen);
    histDecLengthHits->SetLineWidth(2);


    c->cd();
    histDecLengthRec->Draw();
    histDecLengthRecCorrect->Draw("same");
    histDecLengthHits->Draw("same");
    auto legend = new TLegend(0.55, 0.2, 0.85, 0.4);
    legend->SetMargin(0.10);
    legend->SetTextSize(0.03);

    legend->AddEntry(histDecLengthRec, "Generated decay length for ITS reco tracks");
    legend->AddEntry(histDecLengthRecCorrect, "Generated decay length for ITS reco tracks (w/o fake hits)");
    legend->AddEntry(histDecLengthHits, "Generated decay length for tracks with ITS hits");
    legend->Draw();

    c->Write();
    histFakeHits->Write();
    outFile.Close();
};

double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
{
    auto idStart = motherTrack.getFirstDaughterTrackId();
    auto idStop = motherTrack.getLastDaughterTrackId();
    for (auto iD{idStart}; iD < idStop; ++iD)
    {
        auto dauTrack = MCTracks->at(iD);
        if (dauTrack.GetPdgCode() == dauPDG)
        {
            auto decLength = TMath::Sqrt((dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) * (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) + (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) * (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()));
            return decLength;
        }
    }
    return -1;
}