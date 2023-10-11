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

    # four-vector sum of first four elements of each
    # lepton collection (possibly fewer)
    dielectron = events.Electron[:, :4].sum(axis=1)
    dimuon = events.Muon[:, :4].sum(axis=1)

    # sum the results to form the four-lepton four-vector
    fourlep = dielectron + dimuon

    # total number of leptons per event
    n_leptons = (
        ak.num(events.Electron, axis=1) +
        ak.num(events.Muon, axis=1)
    )

    # four-lepton mass, taking into account only events with at least four leptons,
    # and otherwise substituting a predefined EMPTY_FLOAT value
    fourlep_mass = ak.where(
        n_leptons >= 4,
        fourlep.mass,
        EMPTY_FLOAT,
    )

    # write out the resulting mass to the `events` array,
    events = set_ak_column_f32(
        events,
        "m4l",
        fourlep_mass,
    )

    # return the events
    return events
