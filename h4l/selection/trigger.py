
from __future__ import annotations

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import

np = maybe_import("numpy")
ak = maybe_import("awkward")


@selector
def trigger_selection(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:

    # start with an all-false mask
    sel_trigger = ak.Array(np.zeros(len(events), dtype=bool))

    # pick events that passed one of the required triggers
    for trigger in self.dataset_inst.x("require_triggers", []):
        sel_trigger = sel_trigger | events.HLT[trigger]

    # but reject events that also passed one of the triggers to veto
    for trigger in self.dataset_inst.x("veto_triggers", []):
        sel_trigger = sel_trigger & ~events.HLT[trigger]

    return events, SelectionResult(
        steps={
            "trigger": sel_trigger,
        },
    )


@trigger_selection.init
def trigger_selection_init(self: Selector) -> None:
    # return immediately if config object has not been loaded yet
    if not getattr(self, "config_inst", None):
        return

    # add HLT trigger bits to uses
    self.uses |= {
        f"HLT.{trigger}"
        for trigger in self.config_inst.x.all_triggers
    }
