#if !defined(CLING) || defined(ROOTCLING)
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "CommonDataFormat/RangeReference.h"
#include "DetectorsVertexing/PVertexerHelpers.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITS/TrackITS.h"
#include "ReconstructionDataFormats/DCA.h"
#include "TMath.h"
#endif

using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;
using GIndex = o2::dataformats::VtxTrackIndex;
using RRef = o2::dataformats::RangeReference<int, int>;

void vert()
{

	auto fVert = TFile::Open("o2_primary_vertex.root");
	auto treeVert = (TTree *)fVert->Get("o2sim");

	auto fTracks = TFile::Open("o2trac_its.root");
	auto treeTrack = (TTree *)fTracks->Get("o2sim");

	std::vector<o2::dataformats::PrimaryVertex> *verts = nullptr;
	std::vector<Vertex> *seedingVerts = nullptr;

	treeVert->SetBranchAddress("PrimaryVertex", &verts);
	treeTrack->SetBranchAddress("Vertices", &seedingVerts);

	auto histSeedVertsZ = new TH1F("histSeedingVert", "Seeding primary vertices Z;Z (cm);counts", 300, -100, 100);
	auto histVertsZ = new TH1F("histVert", "Primary vertices Z;Z (cm);counts", 300, -100, 100);
	auto histVertsXY = new TH2F("histVertXY", "Primary vertices XY;X (cm);Y (cm)", 300, -2, 2, 300, -3, 3);
	auto histXvsZ = new TH2F("histXvsZ", "Primary vertices XZ correlation;Z (cm);X (cm)", 300, -80, 80, 300, -3, 3);
	auto histYvsZ = new TH2F("histYvsZ", "Primary vertices YZ correlation;Z (cm);Y (cm)", 300, -80, 80, 300, -3, 3);

	histVertsZ->SetFillColor(kBlue);
	histVertsZ->SetFillStyle(3003);

	for (auto ent{0}; ent < treeVert->GetEntriesFast(); ++ent)
	{
		treeVert->GetEntry(ent);
		treeTrack->GetEntry(ent);
		for (auto &v : *verts)
		{
			histVertsZ->Fill(v.getZ());
			histVertsXY->Fill(v.getX(), v.getY());
			histXvsZ->Fill(v.getZ(), v.getX());
			histYvsZ->Fill(v.getZ(), v.getY());
		}
		for (auto &sv : *seedingVerts)
		{
			histSeedVertsZ->Fill(sv.getZ());
		}
	}

	gStyle->SetOptFit(11);
	auto fit = histVertsZ->Fit("gaus", "L", "", 15, 65);

	TCanvas *c1 = new TCanvas("vertices", "vertices", 2200, 1000);
	c1->Divide(2, 3);
	c1->cd(4);
	histVertsZ->Draw();
	c1->cd(2);
	histVertsXY->Draw("colz");
	c1->cd(3);
	histXvsZ->Draw("colz");
	c1->cd(1);
	histYvsZ->Draw("colz");
	c1->cd(5);
	histSeedVertsZ->Draw();
	c1->SaveAs("vertices.png");
	fVert->Close();
}

void dca()
{
	TH2F *histDCApt = new TH2F("dcaPT", "DCA vs #it{p}_{T};#it{p}_{T};DCA (cm)", 150, 0, 10, 150, 0, 2);
	auto fVert = TFile::Open("o2_primary_vertex.root");
	auto fTracks = TFile::Open("o2trac_its.root");

	auto treeVert = (TTree *)fVert->Get("o2sim");
	auto treeTrack = (TTree *)fTracks->Get("o2sim");

	std::vector<o2::dataformats::PrimaryVertex> *verts = nullptr;
	std::vector<o2::vertexing::PVertex> *pvArr = nullptr;
	std::vector<o2::vertexing::V2TRef> *pvRefs = nullptr;
	std::vector<o2::dataformats::VtxTrackIndex> *pvIdx = nullptr;

	treeVert->SetBranchAddress("PrimaryVertex", &verts);
	treeVert->SetBranchAddress("PrimaryVertex", &pvArr);
	treeVert->SetBranchAddress("PV2TrackRefs", &pvRefs);
	treeVert->SetBranchAddress("PVTrackIndices", &pvIdx);

	std::vector<o2::itsmft::ROFRecord> *rofArr = nullptr;
	std::vector<o2::its::TrackITS> *tracks = nullptr;

	treeTrack->SetBranchAddress("ITSTrack", &tracks);
	treeTrack->SetBranchAddress("ITSTracksROF", &rofArr);

	std::vector<std::vector<Vertex>> pvROFs(rofArr->size());
	for (int frame = 0; frame < treeTrack->GetEntriesFast(); frame++)
	{
		if (!treeTrack->GetEvent(frame) || !treeVert->GetEvent(frame))
			continue;
		if ((frame % 100) == 0)
		{
			std::cout << "Processing TF " << frame << std::endl;
		}
		std::vector<std::vector<o2::dataformats::PrimaryVertex>> pvROFs(rofArr->size());
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
					pvROFs[rofId].push_back(pv);
				}
			}
		}

		for (size_t rofId{0}; rofId < rofArr->size(); ++rofId)
		{
			auto &rof = rofArr->at(rofId);

			int trackStart = rof.getFirstEntry();
			int trackStop = trackStart + rof.getNEntries();
			int vcount{0};
			for (auto &v : pvROFs[rofId])
			{
				for (int i = trackStart; i < trackStop; i++)
				{
					auto &track = (*tracks)[i];

					o2::dataformats::DCA dca;
					track.propagateToDCA(v, 2, &dca);
					histDCApt->Fill(track.getPt(), TMath::Sqrt(dca.getR2()));
				}
			}
		}
	}
	auto c1 = new TCanvas("dca", "dca", 1200, 800);
	c1->cd();
	histDCApt->Draw("colz");
	c1->SaveAs("dca.png");
}

void plots()
{
	// gStyle->SetOptStat(0);
	vert();
	dca();
}
