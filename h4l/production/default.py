# coding: utf-8

import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, dev_sandbox
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.production.util import attach_coffea_behavior
from columnflow.production.cms.electron import electron_weights
from columnflow.production.cms.muon import muon_weights
from columnflow.production.categories import category_ids
from columnflow.production.normalization import normalization_weights

from h4l.production.invariant_mass import four_lep_invariant_mass
from h4l.production.weights import (normalized_murmuf_weight,
                                    # normalized_pdf_weight,
                                    normalized_pu_weight)

np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

@producer(
    uses={attach_coffea_behavior, four_lep_invariant_mass, normalized_pu_weight,
          normalized_murmuf_weight, electron_weights, category_ids,
          muon_weights, normalization_weights,
          },
    produces={attach_coffea_behavior, four_lep_invariant_mass, normalized_pu_weight,
          normalized_murmuf_weight, electron_weights, category_ids,
          muon_weights, normalization_weights,
          },
    sandbox=dev_sandbox("bash::$CF_BASE/sandboxes/venv_columnar.sh"),

)
def default(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    events = self[category_ids](events, **kwargs)

    if self.dataset_inst.is_mc:
        # normalization weights
        events = self[normalization_weights](events, **kwargs)

        # normalized pdf weight
        # events = self[normalized_pdf_weight](events, **kwargs)

        # normalized renorm./fact. weight
        events = self[normalized_murmuf_weight](events, **kwargs)

        # normalized pu weights
        events = self[normalized_pu_weight](events, **kwargs)

        # electron weights
        events = self[electron_weights](events, **kwargs)

        # muon weights
        events = self[muon_weights](events, **kwargs)

    events = self[four_lep_invariant_mass](events, **kwargs)

    return events