
"""
Calibration methods.
"""

from columnflow.calibration import Calibrator, calibrator
from columnflow.util import maybe_import

from h4l.calibration.jets import jet_energy

ak = maybe_import("awkward")


@calibrator(
    uses={jet_energy},
    produces={jet_energy},
)
def default(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    """Default calibrator. Corrects only JECs for now."""

    events = self[jet_energy](events, **kwargs)

    return events
