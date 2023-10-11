# coding: utf-8

"""
Jet energy calibration methods.
"""

from columnflow.calibration import Calibrator, calibrator
from columnflow.calibration.cms.jets import jec, jer
from columnflow.util import maybe_import

ak = maybe_import("awkward")


# custom jec calibrator that only runs nominal correction
jec_nominal = jec.derive("jec_nominal", cls_dict={"uncertainty_sources": []})


@calibrator(
    uses={jec_nominal},
    produces={jec_nominal},
)
def jet_energy(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    """
    Common calibrator for jet energy corrections, applying nominal JEC for data, and JEC with
    uncertainties plus JER for MC. Information about used and produced columns and dependent
    calibrators is added in a custom init function below.
    """
    # correct jet energy scale
    events = self[jec_nominal](events, **kwargs)

    # jet energy resolution smearing (MC only)
    if self.dataset_inst.is_mc:
        events = self[jer](events, **kwargs)

    return events


@jet_energy.init
def jet_energy_init(self: Calibrator) -> None:
    # return immediately if dataset object has not been loaded yet
    if not getattr(self, "dataset_inst", None):
        return

    # add columns producs by JER smearing calibrator (MC only)
    if self.dataset_inst.is_mc:
        self.uses.add(jer)
        self.produces.add(jer)
