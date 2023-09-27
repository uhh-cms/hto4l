# coding: utf-8

import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.production.util import attach_coffea_behavior

np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


@producer(
    uses=(
        # "nJet", "Jet.pt", "Jet.eta", "Jet.phi", "Jet.mass", "Jet.area", 
        # "Jet.rawFactor", "Jet.jetId",
        
        {f"{field}.{var}"
         for field in ["Electron", "Muon"]
         for var in ["pt", "mass", "eta", "phi", "charge"]
         }|
         {
             attach_coffea_behavior,
         }
    ),
    produces={
        "m4l",
    },
)
def four_lep_invariant_mass(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    # category ids
    # from IPython import embed
    # embed()
    events = self[attach_coffea_behavior](
        events,
        collections=["Electron", "Muon"],
        **kwargs
    )
    from IPython import embed; embed()

    diele = events.Electron[:, :2].sum(axis=1)
    dimuon = events.Muon[:, :2].sum(axis=1)
    fourlep = diele + dimuon
    diele_mask = ak.num(events.Electron, axis=1) >= 2
    dimuon_mask = ak.num(events.Muon, axis=1) >= 2
    events = set_ak_column_f32(
        events,
        "m4l",
        ak.where(diele_mask&dimuon_mask, fourlep.mass, EMPTY_FLOAT),
    )

    return events
