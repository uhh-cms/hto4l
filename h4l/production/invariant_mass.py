# coding: utf-8

import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, set_ak_column
from columnflow.production.util import attach_coffea_behavior

np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


@producer(
    uses=(
        # "nJet", "Jet.pt", "Jet.eta", "Jet.phi", "Jet.mass", "Jet.area",
        # "Jet.rawFactor", "Jet.jetId",
        {
            f"{field}.{var}"
            for field in ["Electron", "Muon"]
            for var in ["pt", "mass", "eta", "phi", "charge"]
        } | {
            attach_coffea_behavior,
        }
    ),
    produces={
        "m4l",
    },
)
def four_lep_invariant_mass(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Construct four-lepton invariant mass given the Electron and Muon arrays.
    """

    # attach coffea behavior for four-vector arithmetic
    events = self[attach_coffea_behavior](
        events,
        collections=["Electron", "Muon"],
        **kwargs,
    )

    # pad the lepton arrays with `None` entries
    # to ensure they have at least length 2
    electron = ak.pad_none(events.Electron, 2, axis=1)
    muon = ak.pad_none(events.Muon, 2, axis=1)

    # sum over first two elements of each collection
    # (`None` entries propagate)
    dielectron = electron[:, :2].sum(axis=1)
    dimuon = muon[:, :2].sum(axis=1)

    # sum the results to form the four-lepton four-vector
    fourlep = dielectron + dimuon

    # write out the resulting mass to the `events` array,
    # replacing `None` values with a predefined EMPTY_FLOAT value
    events = set_ak_column_f32(
        events,
        "m4l",
        ak.fill_none(fourlep.mass, EMPTY_FLOAT),
    )

    # return the events
    return events
