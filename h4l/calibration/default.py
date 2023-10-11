"""
Default corrections for h4l analysis
"""

from columnflow.calibration import Calibrator, calibrator
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column
from columnflow.production.cms.mc_weight import mc_weight


from h4l.calibration.fsr_photon_correction import fsr_photon_calibrator

np = maybe_import("numpy")
ak = maybe_import("awkward")


@calibrator(
    uses={fsr_photon_calibrator, mc_weight, deterministic_seeds},
    produces={fsr_photon_calibrator, mc_weight, deterministic_seeds},
    exposed=True,
)
def default(
    self: Calibrator,
    events: ak.Array,
    **kwargs,
):
    if self.dataset_inst.is_mc:
        events = self[mc_weight](events, **kwargs)

    events = self[deterministic_seeds](events, **kwargs)

    events = self[fsr_photon_calibrator](events, **kwargs)
    return events
