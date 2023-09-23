
from __future__ import annotations

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import set_ak_column
from columnflow.util import DotDict, maybe_import, dev_sandbox
from h4l.util import IF_NANO_V9, IF_NANO_V10

np = maybe_import("numpy")
ak = maybe_import("awkward")

@selector(
    uses={f"Electron.{var}"
          for var in {"pt", "eta", "deltaEtaSC",
                      "dxy", "dz", "sip3d"}
          } | {
            IF_NANO_V9("Electron.mvaFall17V2Iso"),
            IF_NANO_V10("Electron.mvaHZZIso"),
          },
    exposed=False,
    sandbox=dev_sandbox("bash::$CF_BASE/sandboxes/venv_columnar.sh"),
)
def electron_selection(
    self: Selector,
    events: ak.Array,
    **kwargs,
):
    min_pt = 15
    pt = events.Electron.pt
    fSCeta = abs(events.Electron.eta + events.Electron.deltaEtaSC)

    def return_cuts():

        era = self.config_inst.campaign.x.year
        if era == 2017 or era == 2018 :
            if self.config_inst.campaign.x.version < 10:
                BDT = events.Electron.mvaFall17V2Iso     # Using 2017 WP and training (ElectronMVAEstimatorRun2Fall17IsoV2Values) since this is the only one available in Run2 UL nanoAODs.
            else :
                BDT = events.Electron.mvaHZZIso
            if self.config_inst.campaign.x.preUL : # pre-UL WP for Run II (miniAOD branch: Run2_CutBased_BTag16)
                # print("This is preUL!")
                # from IPython import embed; embed()
                return ((pt<=10.) &     ((fSCeta<0.8                   & BDT > 0.85216885148) | \
                                        (fSCeta>=0.8 & fSCeta<1.479 & BDT > 0.82684550976) | \
                                        (fSCeta>=1.479                & BDT > 0.86937630022))) \
                        | ((pt>10.) &  ((fSCeta<0.8                   & BDT > 0.98248928759) | \
                                        (fSCeta>=0.8 & fSCeta<1.479 & BDT > 0.96919224579) | \
                                        (fSCeta>=1.479                & BDT > 0.79349796445)))

            else: # UL WP (miniAOD branch Run2_CutBased_UL)
                # print("This is not preUL!")
                # from IPython import embed; embed()
                return ((pt<=10.) &     (((fSCeta<0.8)                   & (BDT > 0.9128577458)) | \
                                        ((fSCeta>=0.8) & (fSCeta<1.479) & (BDT > 0.9056792368)) | \
                                        ((fSCeta>=1.479)                & (BDT > 0.9439440575)))) \
                        | ((pt>10.) &  (((fSCeta<0.8)                   & (BDT > 0.1559788054)) | \
                                        ((fSCeta>=0.8) & (fSCeta<1.479) & (BDT > 0.0273863727)) | \
                                        ((fSCeta>=1.479)                & (BDT > -0.5532483665))))
        return ak.ones_like(events.Electron.pt)
    default_mask = (
        (pt > min_pt) &
        (abs(events.Electron.eta) < 2.5) & 
        (events.Electron.dxy < 0.5) &
        (events.Electron.dz < 1.0) & 
        (abs(events.Electron.sip3d) < 4) & True
        # return_cuts()
    )
    electron_idx = ak.local_index(events.Electron.pt, axis=1)
    selected_electron_idx = electron_idx[default_mask]
    # sort for pt
    selected_electron_idx = ak.argsort(
        events.Electron[selected_electron_idx].pt,
        axis=1,
        ascending=False
    )
    


    return events, SelectionResult(
        objects = {
            "Electron": {
                "Electron": selected_electron_idx
            }
        }
    )

@selector(
    uses={f"Muon.{var}" for var in ["pt", "eta", "phi", "mass"]},
    exposed=False,
)
def muon_selector(self: Selector, events: ak.Array, **kwargs) -> tuple[ak.Array, SelectionResult]:
    min_pt = 20
    max_eta = 2.4
    sorted_muons = ak.argsort(events.Muon.pt, axis=1, ascending=False)
    selected_muon_mask = (
        (events.Muon.pt[sorted_muons] > min_pt) & 
        (abs(events.Muon.eta[sorted_muons]) < max_eta)
    )
    muon_idx = sorted_muons[selected_muon_mask]
    return events, SelectionResult(
        objects= {
            "Muon": {
                "Muon": muon_idx
            }
        }
    )