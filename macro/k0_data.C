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

void k0_data(std::string secFile = "o2trac_its.root", std::string pvFileName = "o2_primary_vertex.root")
{

	auto fSecondaries = TFile::Open(secFile.data());
	auto recTree = (TTree *)fSecondaries->Get("o2sim");

	std::vector<o2::its::TrackITS> *tracks = nullptr;

	recTree->SetBranchAddress("ITSTrack", &tracks);
	std::vector<o2::itsmft::ROFRecord> *rofArr = nullptr;
	recTree->SetBranchAddress("ITSTracksROF", &rofArr);

	auto pvFile = TFile::Open(pvFileName.data());
	auto pvTree = (TTree*)pvFile->Get("o2sim");

	std::vector<o2::vertexing::PVertex> *pvArr{nullptr};
	std::vector<o2::vertexing::V2TRef> *pvRefs{nullptr};
	std::vector<o2::dataformats::VtxTrackIndex> *pvIdx{nullptr};


	pvTree->SetBranchAddress("PrimaryVertex", &pvArr);
	pvTree->SetBranchAddress("PV2TrackRefs", &pvRefs);
	pvTree->SetBranchAddress("PVTrackIndices", &pvIdx);

	auto outFile = TFile("k0_inv_mass.root", "recreate");
	TH1D *histInvMass = new TH1D("k0", "k0s rec mass", 80, 0.4, 0.6);
	TH1D *histInvMassLS = new TH1D("LS", "rec LS", 80, 0.4, 0.6);
	TH1D *cosPAH = new TH1D("cosPA", ";cos(#theta_{P})", 1000, 0.9, 1.);

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
		for (size_t idx{0}; idx < pvArr->size(); ++idx) {
			auto& pv{pvArr->at(idx)};
			auto& refs{pvRefs->at(idx)};
			int first{refs.getFirstEntryOfSource(0)}; // only ITS business here
			int last{first + refs.getEntriesOfSource(0) - 1};
			for (size_t rofId{0}; rofId < rofArr->size(); ++rofId)
			{
				auto &rof = rofArr->at(rofId);

				int trackStart = rof.getFirstEntry();
				int trackStop = trackStart + rof.getNEntries();
				int count{0};
				for (int i{first}; i < last; ++i) {
					int current{pvIdx->at(i)};
					if (current >= trackStop || current < trackStart) {
						count--;	
					} else {
						count++;
					}
				}
				if (count > 0) {
					pvROFs[rofId].push_back({pv.getX(), pv.getY(), pv.getZ()});
				}

			}
		}

		for (size_t rofId{0}; rofId < rofArr->size(); ++rofId)
		{
			auto &rof = rofArr->at(rofId);

			int trackStart = rof.getFirstEntry();
			int trackStop = trackStart + rof.getNEntries();

			for (int i = trackStart; i < trackStop; i++)
			{
				auto &track1 = (*tracks)[i];
				for (int j = i + 1; j < trackStop; j++)
				{
					auto &track2 = (*tracks)[j];

					int nCand = mFitterV0.process(track1, track2);
					if (!nCand)
						continue;
					mFitterV0.propagateTracksToVertex();
					std::array<float, 3> p;
					TLorentzVector moth, prong;

					Vec3 pMom{0.f,0.f,0.f};
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

					auto sv = mFitterV0.getPCACandidate();
					bool survives{false};
					for (auto& pv : pvROFs[rofId]) {
						auto delta = pv - sv;
						double cosPA{ROOT::Math::Dot(delta,pMom) / ROOT::Math::Mag(delta) / ROOT::Math::Mag(pMom)};
						cosPAH->Fill(cosPA);
						survives = cosPA > 0.999;
						if (survives) break;
					}
					if (!survives) continue;


					if (track1.getSign() != track2.getSign())
						histInvMass->Fill(moth.M());
					else
						histInvMassLS->Fill(moth.M());
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
	std::cout<<integral/integralLS<<std::endl;
	histInvMassLS->Scale(integral / integralLS);
	histInvMassLS->Draw("same");
	cv.Write();
	cosPAH->Write();
	outFile.Close();
}