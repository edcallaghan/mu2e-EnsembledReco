// Ed Callaghan
// Temporary ntuple containing ensembled collections of TrkInfos
// July 2025

// art
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// art_root_io
#include "art_root_io/TFileService.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// root
#include "TTree.h"

// mu2e
#include "EnsembledReco/inc/KalSeedCollectionCollection.hh"
#include "EventNtuple/inc/InfoStructHelper.hh"
#include "EventNtuple/inc/TrkInfo.hh"
#include "EventNtuple/inc/TrkSegInfo.hh"

namespace mu2e{
  class EnsembledNtupler: public art::EDAnalyzer{
    public:
      struct Config{
        fhicl::Atom<art::InputTag> tag{
          fhicl::Name("KalSeedCollectionCollection"),
          fhicl::Comment("KalSeedCollectionCollection to ntuple")
        };
      };

      using Parameters = art::EDAnalyzer::Table<Config>;
      EnsembledNtupler(const Parameters&);

      virtual void beginJob();
      virtual void analyze(const art::Event&);

    protected:
      art::InputTag _tag;
      std::vector<std::vector<std::vector<TrkInfo>>> _trkinfos;
      std::vector<std::vector<std::vector<std::vector<TrkSegInfo>>>> _seginfos;
      TTree* _tree;

    private:
      /**/
  };

  EnsembledNtupler::EnsembledNtupler(const Parameters& config):
      art::EDAnalyzer(config),
      _tag(config().tag()){
    // framework hooks
    this->consumes<KalSeedCollectionCollection>(_tag);
  }

  void EnsembledNtupler::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    _tree = tfs->make<TTree>("TPerTrack", "TPerTrack");
    _tree->Branch("TrkInfo", &_trkinfos, 32000, 99);
    _tree->Branch("TrkSegInfo", &_seginfos);
  }

  void EnsembledNtupler::analyze(const art::Event& event){
    InfoStructHelper helper;
    const auto handle = event.getValidHandle<KalSeedCollectionCollection>(_tag);

    _trkinfos.clear();
    _seginfos.clear();
    for (const auto& collection: (*handle)){
      std::vector<std::vector<TrkInfo>> tmptrks;
      std::vector<std::vector<std::vector<TrkSegInfo>>> tmpsegs;
      for (const auto& kalseed: collection){
        std::vector<TrkInfo> trks;
        helper.fillTrkInfo(kalseed, trks);
        tmptrks.push_back(trks);

        std::vector<std::vector<TrkSegInfo>> segs;
        helper.fillTrkSegInfo(kalseed, segs);
        tmpsegs.push_back(segs);
      }
        _trkinfos.push_back(tmptrks);
        _seginfos.push_back(tmpsegs);
    }

    _tree->Fill();
  };
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EnsembledNtupler)
