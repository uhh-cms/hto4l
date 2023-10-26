import os
import itertools
import functools

import yaml
import law
import order as od
from scinum import Number

from columnflow.util import DotDict, dev_sandbox
from columnflow.config_util import (
    get_root_processes_from_campaign, add_shift_aliases, get_shifts_from_sources,
    verify_config_processes, add_category,
)


def add_podas_config(
    analysis: od.Analysis,
    campaign: od.Campaign,
    config_name: str or None = None,
    config_id: int or None = None,
    limit_dataset_files: int or None = None,
) -> od.Config:

    # get all root processes
    procs = get_root_processes_from_campaign(campaign)

    # create a config by passing the campaign, so id and name will be identical
    cfg = analysis.add_config(campaign, name=config_name, id=config_id)

    # gather campaign data
    year = campaign.x.year

    corr_postfix = f"{campaign.x.vfp}VFP" if year == 2016 else ""

    # add processes we are interested in
    process_names = [
        # "data_mu",
        "h_ggf_4l",
    ]
    for process_name in process_names:
        # add the process
        proc = cfg.add_process(procs.get(process_name))

        # configuration of colors, labels, etc. can happen here
        if proc.is_mc:
            proc.color1 = (244, 182, 66) if proc.name == "tt" else (244, 93, 66)

    # add datasets we need to study
    dataset_names = [
        # data
        # "data_mu_a",
        # backgrounds
        # signals
        "h_ggf_4l_powheg",
    ]
    for dataset_name in dataset_names:
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))

        # for testing purposes, limit the number of files
        if isinstance(limit_dataset_files, int) and limit_dataset_files > 0:
            for info in dataset.info.values():
                info.n_files = min(info.n_files, limit_dataset_files)

    # verify that the root process of all datasets is part of any of the registered processes
    verify_config_processes(cfg, warn=True)

    # default objects, such as calibrator, selector, producer, ml model, inference model, etc
    cfg.x.default_calibrator = "default"
    cfg.x.default_selector = "default"
    cfg.x.default_producer = "default"
    cfg.x.default_ml_model = None
    cfg.x.default_inference_model = "example"
    cfg.x.default_categories = ("incl",)
    cfg.x.default_variables = ("n_jet", "jet1_pt")

    # process groups for conveniently looping over certain processs
    # (used in wrapper_factory and during plotting)
    cfg.x.process_groups = {}

    # dataset groups for conveniently looping over certain datasets
    # (used in wrapper_factory and during plotting)
    cfg.x.dataset_groups = {}

    # category groups for conveniently looping over certain categories
    # (used during plotting)
    cfg.x.category_groups = {}

    # variable groups for conveniently looping over certain variables
    # (used during plotting)
    cfg.x.variable_groups = {}

    # shift groups for conveniently looping over certain shifts
    # (used during plotting)
    cfg.x.shift_groups = {}

    # selector step groups for conveniently looping over certain steps
    # (used in cutflow tasks)
    cfg.x.selector_step_groups = {
        "default": ["muon", "jet"],
    }

    # custom method and sandbox for determining dataset lfns
    cfg.x.get_dataset_lfns = None
    cfg.x.get_dataset_lfns_sandbox = None

    # whether to validate the number of obtained LFNs in GetDatasetLFNs
    # (currently set to false because the number of files per dataset is truncated to 2)
    cfg.x.validate_dataset_lfns = False

    # lumi values in inverse pb
    # https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2?rev=2#Combination_and_correlations
    cfg.x.luminosity = Number(41480, {
        "lumi_13TeV_2017": 0.02j,
        "lumi_13TeV_1718": 0.006j,
        "lumi_13TeV_correlated": 0.009j,
    })

    # names of muon correction sets and working points
    # (used in the muon producer)
    cfg.x.muon_sf_names = ("NUM_TightRelIso_DEN_TightIDandIPCut", f"{year}_UL")

    # register shifts
    cfg.add_shift(name="nominal", id=0)

    # tune shifts are covered by dedicated, varied datasets, so tag the shift as "disjoint_from_nominal"
    # (this is currently used to decide whether ML evaluations are done on the full shifted dataset)
    cfg.add_shift(name="tune_up", id=1, type="shape", tags={"disjoint_from_nominal"})
    cfg.add_shift(name="tune_down", id=2, type="shape", tags={"disjoint_from_nominal"})

    # fake jet energy correction shift, with aliases flaged as "selection_dependent", i.e. the aliases
    # affect columns that might change the output of the event selection
    cfg.add_shift(name="jec_up", id=20, type="shape")
    cfg.add_shift(name="jec_down", id=21, type="shape")
    add_shift_aliases(
        cfg,
        "jec",
        {
            "Jet.pt": "Jet.pt_{name}",
            "Jet.mass": "Jet.mass_{name}",
            "MET.pt": "MET.pt_{name}",
            "MET.phi": "MET.phi_{name}",
        },
    )

    # event weights due to muon scale factors
    cfg.add_shift(name="mu_up", id=10, type="shape")
    cfg.add_shift(name="mu_down", id=11, type="shape")
    add_shift_aliases(cfg, "mu", {"muon_weight": "muon_weight_{direction}"})

    cfg.add_shift(name="ele_up", id=12, type="shape")
    cfg.add_shift(name="ele_down", id=13, type="shape")
    add_shift_aliases(cfg, "ele", {"electron_weight": "electron_weight_{direction}"})

    cfg.add_shift(name="murmuf_up", id=140, type="shape")
    cfg.add_shift(name="murmuf_down", id=141, type="shape")
    add_shift_aliases(
        cfg,
        "murmuf",
        {
            "murmuf_weight": "murmuf_weight_{direction}",
            "normalized_murmuf_weight": "normalized_murmuf_weight_{direction}",
        },
    )

    # external files
    json_mirror = "/afs/cern.ch/user/m/mrieger/public/mirrors/jsonpog-integration-dfd90038"
    cfg.x.external_files = DotDict.wrap({
        # jet energy correction
        "jet_jerc": (f"{json_mirror}/POG/JME/{year}{corr_postfix}_UL/jet_jerc.json.gz", "v1"),

        # tau energy correction and scale factors
        "tau_sf": (f"{json_mirror}/POG/TAU/{year}{corr_postfix}_UL/tau.json.gz", "v1"),

        # electron scale factors
        "electron_sf": (f"{json_mirror}/POG/EGM/{year}{corr_postfix}_UL/electron.json.gz", "v1"),

        # muon scale factors
        "muon_sf": (f"{json_mirror}/POG/MUO/{year}{corr_postfix}_UL/muon_Z.json.gz", "v1"),

        # btag scale factor
        "btag_sf_corr": (f"{json_mirror}/POG/BTV/{year}{corr_postfix}_UL/btagging.json.gz", "v1"),

        # met phi corrector
        "met_phi_corr": (f"{json_mirror}/POG/JME/{year}{corr_postfix}_UL/met.json.gz", "v1"),

        # hh-btag repository (lightweight) with TF saved model directories
        "hh_btag_repo": ("https://github.com/hh-italian-group/HHbtag/archive/1dc426053418e1cab2aec021802faf31ddf3c5cd.tar.gz", "v1"),  # noqa
    })

    # external files with more complex year dependence
    if year == 2016:
        cfg.x.external_files.update(DotDict.wrap({
            # lumi files
            "lumi": {
                "golden": ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt", "v1"),  # noqa
                "normtag": ("/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json", "v1"),
            },

            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData?rev=45#Pileup_JSON_Files_For_Run_II
            "pu": {
                "json": ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt", "v1"),  # noqa
                "mc_profile": ("https://raw.githubusercontent.com/cms-sw/cmssw/a65c2e1a23f2e7fe036237e2e34cda8af06b8182/SimGeneral/MixingModule/python/mix_2016_25ns_UltraLegacy_PoissonOOTPU_cfi.py", "v1"),  # noqa
                "data_profile": {
                    "nominal": (f"/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2016-{campaign.x.vfp}VFP-69200ub-99bins.root", "v1"),  # noqa
                    "minbias_xs_up": (f"/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2016-{campaign.x.vfp}VFP-72400ub-99bins.root", "v1"),  # noqa
                    "minbias_xs_down": (f"/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2016-{campaign.x.vfp}VFP-66000ub-99bins.root", "v1"),  # noqa
                },
            },
        }))
    elif year == 2017:
        cfg.x.external_files.update(DotDict.wrap({
            # lumi files
            "lumi": {
                "golden": ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt", "v1"),  # noqa
                "normtag": ("/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json", "v1"),
            },

            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData?rev=45#Pileup_JSON_Files_For_Run_II
            "pu": {
                "json": ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/UltraLegacy/pileup_latest.txt", "v1"),  # noqa
                "mc_profile": ("https://raw.githubusercontent.com/cms-sw/cmssw/435f0b04c0e318c1036a6b95eb169181bbbe8344/SimGeneral/MixingModule/python/mix_2017_25ns_UltraLegacy_PoissonOOTPU_cfi.py", "v1"),  # noqa
                "data_profile": {
                    "nominal": ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root", "v1"),  # noqa
                    "minbias_xs_up": ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2017-72400ub-99bins.root", "v1"),  # noqa
                    "minbias_xs_down": ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2017-66000ub-99bins.root", "v1"),  # noqa
                },
            },
        }))
    else:  # year 2018
        cfg.x.external_files.update(DotDict.wrap({
            # lumi files
            "lumi": {
                "golden": ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt", "v1"),  # noqa
                "normtag": ("/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json", "v1"),
            },

            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData?rev=45#Pileup_JSON_Files_For_Run_II
            "pu": {
                "json": ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/pileup_latest.txt", "v1"),  # noqa
                "mc_profile": ("https://raw.githubusercontent.com/cms-sw/cmssw/a65c2e1a23f2e7fe036237e2e34cda8af06b8182/SimGeneral/MixingModule/python/mix_2018_25ns_UltraLegacy_PoissonOOTPU_cfi.py", "v1"),  # noqa
                "data_profile": {
                    "nominal": ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root", "v1"),  # noqa
                    "minbias_xs_up": ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2018-72400ub-99bins.root", "v1"),  # noqa
                    "minbias_xs_down": ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins.root", "v1"),  # noqa
                },
            },
        }))

    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0

    # columns to keep after certain steps
    cfg.x.keep_columns = DotDict.wrap({
        "cf.ReduceEvents": {
            # general event info
            "run", "luminosityBlock", "event",
            # object info
            "Jet.btagDeepFlavB", "Jet.hadronFlavour",
            "Muon.pfRelIso04_all", "Muon.charge",
            "Electron.deltaEtaSC", "Electron.charge",
            "Electron.mvaFall17V2Iso", "Electron.mvaHZZIso",
            "MET.pt", "MET.phi", "MET.significance", "MET.covXX", "MET.covXY",
            "MET.covYY",
            "PV.npvs",
            # columns added during selection
            "deterministic_seed", "process_id", "mc_weight", "cutflow.*",
            "channel_id", "category_ids", "mc_weight", "pdf_weight*", "murmuf_weight*",
            "leptons_os", "single_triggered", "cross_triggered",
            "pu_weight*",
        } | {
            # four momenta information
            f"{field}.{var}"
            for field in ["Jet", "Muon", "Electron"]
            for var in ["pt", "eta", "phi", "mass", "e"]
        } | {
            f"{field}.fsr_uncorrected_{var}"
            for field in ["Muon", "Electron"]
            for var in ["pt", "eta", "phi", "mass", "e"]
        },
        "cf.MergeSelectionMasks": {
            "normalization_weight", "process_id", "category_ids", "cutflow.*",
        },
        "cf.UniteColumns": {
            "*",
        },
    })

    # names of electron correction sets and working points
    # (used in the electron_sf producer)
    cfg.x.electron_sf_names = ("UL-Electron-ID-SF", f"{year}{corr_postfix}", "wp80iso")

    # names of muon correction sets and working points
    # (used in the muon producer)
    cfg.x.muon_sf_names = ("NUM_TightRelIso_DEN_TightIDandIPCut", f"{year}{corr_postfix}_UL")

    # add pt threshold for lepton working points in the H->4l analysis
    # todo: pt for tight WP defined by threshold for scale factors, to be checked
    cfg.x.electron_thresholds = DotDict({"pt": {"loose": 7, "tight": 15}})
    cfg.x.muon_thresholds = DotDict({"pt": {"loose": 5, "tight": 15}})

    # event weight columns as keys in an OrderedDict, mapped to shift instances they depend on
    get_shifts = functools.partial(get_shifts_from_sources, cfg)
    cfg.x.event_weights = DotDict({
        "normalization_weight": [],
        "muon_weight": get_shifts("mu"),
    })

    # versions per task family, either referring to strings or to callables receving the invoking
    # task instance and parameters to be passed to the task family
    cfg.x.versions = {
        # "cf.CalibrateEvents": "prod1",
        # "cf.SelectEvents": (lambda cls, inst, params: "prod1" if params.get("selector") == "default" else "dev1"),
        # ...
    }

    # channels
    # (just one for now)
    cfg.add_channel(name="4mu", id=1)

    # add categories using the "add_category" tool which adds auto-generated ids
    # the "selection" entries refer to names of selectors, e.g. in selection/example.py
    add_category(
        cfg,
        name="incl",
        selection="cat_incl",
        label="inclusive",
    )
    add_category(
        cfg,
        name="2j",
        selection="cat_2j",
        label="2 jets",
    )

    from h4l.config.variables import add_variables
    add_variables(cfg=cfg)