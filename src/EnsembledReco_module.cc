// Ed Callaghan
// Run Kalman Filter reco under various constraint forms, and return weighted sum
// July 2025

// stl
#include <string>
#include <vector>

// art
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// canvas
#include "canvas/Persistency/Common/FindOneP.h"

// cetlib_except
#include "cetlib_except/exception.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "EnsembledReco/inc/KalSeedCollectionCollection.hh"

// maybe stay, maybe go
#include <memory>
#include <vector>
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeneralUtilities/inc/OwningPointerCollection.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/KinKalGeom/inc/SurfaceMap.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateToZ.hh"
#include "Offline/Mu2eKinKal/inc/KKBField.hh"
#include "Offline/Mu2eKinKal/inc/KKConstantBField.hh"
#include "Offline/Mu2eKinKal/inc/KKFit.hh"
#include "Offline/Mu2eKinKal/inc/KKMaterial.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHitCluster.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawXing.hh"
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/KalSeedAssns.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "KinKal/Geometry/Disk.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Trajectory/LoopHelix.hh"

namespace mu2e{
  // TODO Trj -> template<typename T>
  using Prm = KinKal::LoopHelix;
  using Trj = KinKal::ParticleTrajectory<Prm>;
  using Trk = KKTrack<Prm>;
  using KinKalStrawHits = std::vector< std::shared_ptr< KinKal::Hit<Prm> > >;
  using KinKalStrawXngs = std::vector< std::shared_ptr< KinKal::ElementXing<Prm> > >;

  class EnsembledReco: public art::EDProducer{
    public:
      struct Config{
        fhicl::Atom<art::InputTag> tag{
          fhicl::Name("KalSeedCollection"),
          fhicl::Comment("KalSeedCollection to analyze")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      EnsembledReco(const Parameters&);

      void produce(art::Event& event);
      void sample(Trk&, SurfaceMap::SurfacePairCollection&,
                  double, bool, bool) const;
      void sample(Trk&, double, double, double) const;
      void LoopHelixFit_sampleFit(Trk&, SurfaceMap::SurfacePairCollection&,
                                  double, bool, bool) const;
      void LoopHelixFit_toTrackerEnds(Trk&, double, double, double) const;

      using KinKalAuxiliaries = struct {
        KKFit<Prm>& kkfit;
        KKMaterial& kkmat;
        std::shared_ptr<KinKal::BFieldMap> kkbf;
        const ComboHitCollection& hits;
        const StrawResponse& response;
        const Tracker& tracker;
        const Calorimeter& calorimeter;
      };
      using KinKalFittables = struct {
        Prm seed;
        PDGCode::type pdg;
        art::Ptr<HelixSeed> helix;
        std::vector< std::shared_ptr< KKStrawHit  <Prm> > > shits;
        std::vector< std::shared_ptr< KKCaloHit   <Prm> > > chits;
        std::vector< std::shared_ptr< KKStrawXing <Prm> > > sxngs;
      };
      void reconstitute(const Trk&, const KalSeed&,
                        KinKalAuxiliaries&, KinKalFittables&);

    protected:
      art::InputTag _tag;
      ProditionsHandle<StrawResponse> _straw_proditions;
      ProditionsHandle<Tracker> _tracker_proditions;

    private:
      /**/
  };

  EnsembledReco::EnsembledReco(const Parameters& config):
      art::EDProducer(config),
      _tag(config().tag()){
    // framework hooks
    this->consumes<KalSeedCollection>(_tag);
    this->produces<KalSeedCollectionCollection>();
  }

  void EnsembledReco::produce(art::Event& event){
    // for some reason, need a reference to the calorimeter...
    GeomHandle<Calorimeter> calo;

    // fetch seed fits...
    using KKTrackCollection = OwningPointerCollection< KKTrack<Prm> >;
    auto kth = event.getValidHandle<KKTrackCollection>(_tag);
    auto ksh = event.getValidHandle<KalSeedCollection>(_tag);
    if (ksh->size() != kth->size()){
      std::string msg = "KKTrack / KalSeed count mismatch";
      throw cet::exception("EnsembledReco") << msg << std::endl;
    }

    // ...and their configuration
    auto provenance = ksh.provenance();
    auto pset = provenance->parameterSet();

    // forward from config used for original fit
    auto kkfcps = pset.get<fhicl::ParameterSet>("KKFitSettings");
    auto kkfc = fhicl::Table<KKFitConfig>(kkfcps);
    auto kkfit = KKFit<Prm>(kkfc());
    auto fitconfig = Mu2eKinKal::makeConfig(fhicl::Table<Mu2eKinKal::KinKalConfig>(pset.get<fhicl::ParameterSet>("FitSettings"))());

    // for sampling intersections
    auto extconfig = Mu2eKinKal::makeConfig(fhicl::Table<Mu2eKinKal::KinKalConfig>(pset.get<fhicl::ParameterSet>("ExtensionSettings"))());
    auto btol = static_cast<double>(
      pset.get<fhicl::ParameterSet>("ExtensionSettings")
          .get<float>("BCorrTolerance")
    );
    auto maxdt = static_cast<double>(
      pset.get<fhicl::ParameterSet>("Extrapolation").get<float>("MaxDt")
    );

    // for track extensions -.-
    auto modset = pset.get<fhicl::ParameterSet>("ModuleSettings");
    auto cctag = modset.get<art::InputTag>("CaloClusterCollection");
    auto cch = event.getValidHandle<CaloClusterCollection>(cctag);

    // field
    GeomHandle<BFieldManager> bfm;
    GeomHandle<DetectorSystem> dts;
    auto kkbf = std::make_shared<KKBField>(*bfm, *dts);

    // combo hits which were used in the original fit
    auto cht = modset.get<art::InputTag>("ComboHitCollection");
    auto chh = event.getValidHandle<ComboHitCollection>(cht);
    auto hsts = modset.get< std::vector<art::InputTag> >("HelixSeedCollections");
    if (hsts.size() != 1){
      std::string msg = "multiple HelixSeedCollections";
      msg += ", don't know how to proceed.";
      throw cet::exception("EnsembledReco") << msg << std::endl;
    }
    //auto hsc = event.getValidHandle<HelixSeedCollection>(hsts[0]);
    //auto khah = event.getValidHandle<KalHelixAssns>(_tag);
    art::FindOneP<HelixSeed> helices(ksh, event, _tag);

    auto eventid = event.id();
    auto& response = _straw_proditions.get(eventid);
    auto& tracker = _tracker_proditions.get(eventid);

    // downstream configuration for sampling surface intersections
    double intersection_tolerance = modset.get<double>("IntersectionTolerance");
    /*
    bool intersection_in_range = modset.get<bool>("SampleInRange");
    bool intersection_in_bounds = modset.get<bool>("SampleInBounds");
    auto labels = modset.get< std::vector<std::string> >("SampleSurfaces");
    SurfaceIdCollection surfaceids;
    for (const auto& id: labels){
      surfaceids.push_back(SurfaceId(id, -1));
    }
    SurfaceMap::SurfacePairCollection surfaces;
    SurfaceMap().surfaces(surfaceids, surfaces);
    */

    // auxiliaries
    auto kkmat = KKMaterial(fhicl::Table<KKMaterial::Config>(pset.get<fhicl::ParameterSet>("MaterialSettings"))());
    KinKalAuxiliaries auxiliaries{
      kkfit,
      kkmat,
      kkbf,
      //*khah,
      //*hsc,
      *chh,
      response,
      tracker,
      *calo
    };

    // rest come from the KalSeed itself
    KinKalFittables fittables;

    auto out = std::make_unique<KalSeedCollectionCollection>(ksh->size());
    for (size_t i = 0 ; i < ksh->size() ; i++){
      const auto kalseed = ksh->at(i);
      const auto& ktrk = kth->at(i);
      const auto& helix = helices.at(i);
      auto reseed = kkfit.createSeed(ktrk, kalseed.status(), *calo);
      //auto flag = TrkFitFlag();
      auto flag = kalseed.status();

      fittables.helix = helix;
      this->reconstitute(ktrk, kalseed, auxiliaries, fittables);

      // construct track using FitSettings (``seed fit'')
      auto newtrk = KKTrack<Prm>(
        fitconfig,
        *kkbf,
        fittables.seed,
        fittables.pdg,
        auxiliaries.kkfit.strawHitClusterer(),
        fittables.shits,
        fittables.sxngs,
        fittables.chits
      );
      // extend track using ExtensionSettings (``full fit'')
      if (0 < extconfig.schedule().size()){
        auxiliaries.kkfit.extendTrack(
          extconfig,
          *kkbf,
          auxiliaries.tracker,
          auxiliaries.response,
          auxiliaries.kkmat.strawMaterial(),
          auxiliaries.hits,
          auxiliaries.calorimeter,
          cch,
          newtrk
        );
      }
      this->sample(newtrk, intersection_tolerance, maxdt, btol);
      auto newseed = kkfit.createSeed(newtrk, flag, *calo);

      // execute a no-nop fit iteration
//    auto noop = KinKal::Config();
//    noop.schedule_.push_back(KinKal::MetaIterConfig());
      KinKalStrawHits empty_hits;
      KinKalStrawXngs empty_xngs;
//    newtrk.extend(noop, empty_hits, empty_xngs);
//    this->sample(newtrk, intersection_tolerance, maxdt, btol);
//    auto renewseed = kkfit.createSeed(newtrk, flag, *calo);

      // execute a different noop fit iteration
      auto renoop = extconfig;
      renoop.schedule_.clear();
      renoop.schedule_.push_back(KinKal::MetaIterConfig());
      renoop.dwt_ = 1.0;
      newtrk.extend(renoop, empty_hits, empty_xngs);
      this->sample(newtrk, intersection_tolerance, maxdt, btol);
      auto renewseed = kkfit.createSeed(newtrk, flag, *calo);

      (*out)[i].push_back(kalseed);
      (*out)[i].push_back(reseed);
      (*out)[i].push_back(newseed);
      (*out)[i].push_back(renewseed);
    }

    event.put(std::move(out));
  }

  void EnsembledReco::reconstitute(const Trk& track,
                                   const KalSeed& seed,
                                   KinKalAuxiliaries& auxiliaries,
                                   KinKalFittables& fittables){
    // easy ones first
    fittables.seed = *(track.seedTraj().pieces().front());
    fittables.pdg = seed.particle();

    // these require reinterpretation
    std::vector< std::shared_ptr< KKStrawHit<Prm> > > shits;
    std::vector< std::shared_ptr< KKStrawXing<Prm> > > sxngs;
    std::vector< std::shared_ptr< KKCaloHit<Prm> > > chits;

    // find HelixSeed which originally seeded fit

    std::vector<StrawHitIndex> indices;
    fittables.helix->hits().fillStrawHitIndices(indices,
                                                StrawIdMask::uniquestraw);
    shits.reserve(indices.size());
    sxngs.reserve(indices.size());

    auxiliaries.kkfit.makeStrawHits(
      auxiliaries.tracker,
      auxiliaries.response,
      *(auxiliaries.kkbf),
      auxiliaries.kkmat.strawMaterial(),
      fittables.seed,
      auxiliaries.hits,
      indices,
      shits,
      sxngs
    );

    // TODO
    // could fill calo hits here... eh... tbd...
    if (auxiliaries.kkfit.useCalo()){
      auto cluster = fittables.helix->caloCluster();
      if (cluster.isNonnull()){
        auxiliaries.kkfit.makeCaloHit(
          cluster,
          auxiliaries.calorimeter,
          fittables.seed,
          chits
        );
      }
    }

    // store output
    fittables.shits = shits;
    fittables.sxngs = sxngs;
    fittables.chits = chits;
  }

  void EnsembledReco::sample(Trk& track,
                             SurfaceMap::SurfacePairCollection& surfaces,
                             double tolerance,
                             bool in_range,
                             bool in_bounds) const{
    this->LoopHelixFit_sampleFit(track, surfaces,
                                 tolerance, in_range, in_bounds);
  }

  // near-direct copy
  void EnsembledReco::LoopHelixFit_sampleFit(Trk& track,
                              SurfaceMap::SurfacePairCollection& surfaces,
                              double tolerance,
                              bool in_range,
                              bool in_bounds) const{
    // ejc
    auto& kktrk = track;
    auto& sample_ = surfaces;
    auto intertol_ = tolerance;
    auto sampleinrange_ = in_range;
    auto sampleinbounds_ = in_bounds;
    // copy from below
    auto const& ftraj = kktrk.fitTraj();
    std::vector<TimeRange> ranges;
    // test for reflection, and if present, split the test in 2
    auto refltraj = ftraj.reflection(ftraj.range().begin());
    if(refltraj){
      double tmid = refltraj->range().begin();
      ranges.push_back(TimeRange(ftraj.range().begin(),tmid));
      ranges.push_back(TimeRange(tmid,ftraj.range().end()));
    } else {
      ranges.push_back(ftraj.range());
    }
    for(auto range : ranges) {
      double tbeg = range.begin();
      double tend = range.end();
      for(auto const& surf : sample_){
        // search for intersections with each surface within the specified time range, going forwards in time
        bool goodinter(true);
        size_t max_inter = 100; // limit the number of intersections
        size_t cur_inter = 0;
        // loop to find multiple intersections
        while(goodinter && tbeg < tend && cur_inter < max_inter){
          cur_inter += 1;
          TimeRange irange(tbeg,tend);
          auto surfinter = KinKal::intersect(ftraj,*surf.second,irange,intertol_);
          goodinter = surfinter.onsurface_ && ( (! sampleinbounds_) || surfinter.inbounds_ ) && ( (!sampleinrange_) || surfinter.inrange_);
          if(goodinter) {
            // save the intersection information
            kktrk.addIntersection(surf.first,surfinter);
            // update for the next intersection
            // move past existing intersection to avoid repeating
            double epsilon = intertol_/ftraj.speed(surfinter.time_);
            tbeg = surfinter.time_ + epsilon;
          }
        }
      }
    }
  }

  void EnsembledReco::sample(Trk& track,
                             double spatial_tolerance,
                             double maxdt,
                             double field_tolerance) const{
    this->LoopHelixFit_toTrackerEnds(track,
                                     spatial_tolerance,
                                     maxdt,
                                     field_tolerance);
  }

  // near-direct copy
  void EnsembledReco::LoopHelixFit_toTrackerEnds(Trk& track,
                                                 double spatial_tolerance,
                                                 double maxdt,
                                                 double field_tolerance) const{
    // ejc
    auto& ktrk = track;
    auto intertol_ = spatial_tolerance;
    auto trkrep = SurfaceMap().tracker();
    auto trackerFront_ = ExtrapolateToZ(maxdt,
                                        field_tolerance,
                                        trkrep.front().center().Z(),
                                        0);
    auto trackerBack_  = ExtrapolateToZ(maxdt,
                                        field_tolerance,
                                        trkrep.back().center().Z(),
                                        0);
    auto trkfrontptr_ = trkrep.frontPtr();
    auto trkmidptr_ = trkrep.middlePtr();
    auto trkbackptr_ = trkrep.backPtr();
    // copy from below
    // time direction to reach the bounding surfaces from the active region depends on the z momentum. This calculation assumes the particle doesn't
    // reflect inside the tracker volumei
    auto const& ftraj = ktrk.fitTraj();
    auto dir0 = ftraj.direction(ftraj.t0());
    TimeDir fronttdir = (dir0.Z() > 0) ? TimeDir::backwards : TimeDir::forwards;
    TimeDir backtdir = (dir0.Z() > 0) ? TimeDir::forwards : TimeDir::backwards;
    auto tofront = ktrk.extrapolate(fronttdir,trackerFront_);
    auto toback = ktrk.extrapolate(backtdir,trackerBack_);
    // record the standard tracker intersections
    static const SurfaceId tt_front("TT_Front");
    static const SurfaceId tt_mid("TT_Mid");
    static const SurfaceId tt_back("TT_Back");

    // start with the middle
    auto midinter = KinKal::intersect(ftraj,*trkmidptr_,ftraj.range(),intertol_);
    if(midinter.good()) ktrk.addIntersection(tt_mid,midinter);
    if(tofront){
      // check the front piece first; that is usually correct
      // track extrapolation to the front succeeded, but the intersection failed. Use the last trajectory to force an intersection
      auto fhel = fronttdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
      auto frontinter = KinKal::intersect(fhel,*trkfrontptr_,fhel.range(),intertol_,fronttdir);
      if(!frontinter.good()){
        // start from the middle
        TimeRange frange = ftraj.range();
        if(midinter.good())frange = fronttdir == TimeDir::forwards ? TimeRange(midinter.time_,ftraj.range().end()) : TimeRange(ftraj.range().begin(),midinter.time_);
        frontinter = KinKal::intersect(ftraj,*trkfrontptr_,frange,intertol_,fronttdir);
      }
      if(frontinter.good()) ktrk.addIntersection(tt_front,frontinter);
    }
    if(toback){
      // start from the middle
      TimeRange brange = ftraj.range();
      if(midinter.good())brange = backtdir == TimeDir::forwards ? TimeRange(midinter.time_,ftraj.range().end()) : TimeRange(ftraj.range().begin(),midinter.time_);
      auto backinter = KinKal::intersect(ftraj,*trkbackptr_,brange,intertol_,backtdir);
      if(backinter.good())ktrk.addIntersection(tt_back,backinter);
    }
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EnsembledReco)
