# coding: utf-8

"""
Configuration of the hto4l analysis.
"""

import os
import functools

import law
import order as od
from scinum import Number

from columnflow.util import DotDict, maybe_import
from columnflow.columnar_util import EMPTY_FLOAT
from columnflow.config_util import (
    get_root_processes_from_campaign, add_shift_aliases, get_shifts_from_sources, add_category,
    verify_config_processes,
)

ak = maybe_import("awkward")


#
# the main analysis object
#

analysis_h4l = ana = od.Analysis(
    name="analysis_h4l",
    id=1,
)

# analysis-global versions
# (see cfg.x.versions below for more info)
ana.x.versions = {}

# files of bash sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.bash_sandboxes = [
    "$CF_BASE/sandboxes/cf.sh",
    law.config.get("analysis", "default_columnar_sandbox"),
]

# files of cmssw sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.cmssw_sandboxes = [
    # "$CF_BASE/sandboxes/cmssw_default.sh",
]

# clear the list when cmssw bundling is disabled
if not law.util.flag_to_bool(os.getenv("H4L_BUNDLE_CMSSW", "1")):
    del ana.x.cmssw_sandboxes[:]

# config groups for conveniently looping over certain configs
# (used in wrapper_factory)
ana.x.config_groups = {}


#
# setup configs
#

# an example config is setup below, based on cms NanoAOD v9 for Run2 2017, focussing on
# ttbar and single top MCs, plus single muon data
# update this config or add additional ones to accomodate the needs of your analysis

from cmsdb.campaigns.run2_2018_nano_v9 import campaign_run2_2018_nano_v9
from cmsdb.campaigns.run2_2017_nano_v9 import campaign_run2_2017_nano_v9

# copy the campaign
# (creates copies of all linked datasets, processes, etc. to allow for encapsulated customization)
from h4l.config.config_podas import add_podas_config

add_podas_config(
    analysis=analysis_h4l,
    campaign=campaign_run2_2018_nano_v9.copy(),
    config_name=campaign_run2_2018_nano_v9.name,
    config_id=1,
)

# add config for development with limited number of files
add_podas_config(
    analysis=analysis_h4l,
    campaign=campaign_run2_2018_nano_v9.copy(),
    config_name=f"{campaign_run2_2018_nano_v9.name}_limited",
    config_id=2,
    limit_dataset_files=2,
)


add_podas_config(
    analysis=analysis_h4l,
    campaign=campaign_run2_2017_nano_v9.copy(),
    config_name=campaign_run2_2017_nano_v9.name,
    config_id=3,
)

# add config for development with limited number of files
add_podas_config(
    analysis=analysis_h4l,
    campaign=campaign_run2_2017_nano_v9.copy(),
    config_name=f"{campaign_run2_2017_nano_v9.name}_limited",
    config_id=31,
    limit_dataset_files=2,
)