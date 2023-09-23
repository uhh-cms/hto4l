from operator import and_
from functools import reduce
from collections import defaultdict

from columnflow.production.util import attach_coffea_behavior
from columnflow.production.cms.mc_weight import mc_weight
from columnflow.production.cms.pileup import pu_weight
# from columnflow.production.cms.pdf import pdf_weights
from columnflow.production.cms.scale import murmuf_weights
from columnflow.selection.cms.json_filter import json_filter
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.stats import increment_stats
from columnflow.production.processes import process_ids


from columnflow.util import maybe_import, dev_sandbox

from h4l.selection.lepton import electron_selection, muon_selector

np = maybe_import("numpy")
ak = maybe_import("awkward")


@selector(
    uses={"event", attach_coffea_behavior, json_filter, mc_weight, electron_selection,
          pu_weight, murmuf_weights, increment_stats, process_ids, muon_selector,
          },
    produces={attach_coffea_behavior, json_filter, mc_weight, electron_selection,
          pu_weight, murmuf_weights, increment_stats, process_ids, muon_selector,
          },
    sandbox=dev_sandbox("bash::$CF_BASE/sandboxes/venv_columnar.sh"),
    exposed=True,
)
def default(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    # ensure coffea behavior
    events = self[attach_coffea_behavior](events, **kwargs)
    results = SelectionResult()
    results.steps.all = ak.full_like(events.event, True)

    # filter bad data events according to golden lumi mask
    if self.dataset_inst.is_data:
        events, json_filter_results = self[json_filter](events, **kwargs)
        results += json_filter_results

    # add corrected mc weights
    if self.dataset_inst.is_mc:
        events = self[mc_weight](events, **kwargs)

        # mc-only functions
        # pdf weights
        # events = self[pdf_weights](events, **kwargs)

        # renormalization/factorization scale weights
        events = self[murmuf_weights](events, **kwargs)

        # pileup weights
        events = self[pu_weight](events, **kwargs)

    events, ele_results = self[electron_selection](events, call_force=True, **kwargs)
    results += ele_results

    events, muon_results = self[muon_selector](events, call_force=True, **kwargs)
    results += muon_results

    ele_idx = results.objects.Electron.Electron
    muon_idx = results.objects.Muon.Muon
    n_ele = ak.num(events.Electron[ele_idx], axis=1)
    n_muon = ak.num(events.Muon[muon_idx], axis=1)
    results.steps["four_leptons"] = (n_ele + n_muon) >= 4
    event_sel = reduce(and_, results.steps.values())

    events = self[process_ids](events, **kwargs)


    # increment stats
    weight_map = {
        "num_events": Ellipsis,
        "num_events_selected": event_sel,
    }
    group_map = {}
    group_combinations = []
    if self.dataset_inst.is_mc:
        weight_map["sum_mc_weight"] = events.mc_weight
        weight_map["sum_mc_weight_selected"] = (events.mc_weight, event_sel)
        # pu weights with variations
        for name in sorted(self[pu_weight].produces):
            weight_map[f"sum_mc_weight_{name}"] = (events.mc_weight * events[name], Ellipsis)
        # pdf and murmuf weights with variations
        for v in ["", "_up", "_down"]:
            weight_map[f"sum_murmuf_weight{v}"] = events[f"murmuf_weight{v}"]
            weight_map[f"sum_murmuf_weight{v}_selected"] = (events[f"murmuf_weight{v}"], event_sel)
        
        # groups
        group_map = {
            **group_map,
            # per process
            "process": {
                "values": events.process_id,
                "mask_fn": (lambda v: events.process_id == v),
            },
            # # per jet multiplicity
            # "njet": {
            #     "values": results.x.n_central_jets,
            #     "mask_fn": (lambda v: results.x.n_central_jets == v),
            # },
        }
        # combinations
        group_combinations.append(("process",))
    events, results = self[increment_stats](
        events,
        results,
        stats,
        weight_map=weight_map,
        group_map=group_map,
        group_combinations=group_combinations,
        **kwargs,
    )

    return events, results