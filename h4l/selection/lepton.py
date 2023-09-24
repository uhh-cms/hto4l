
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
    sorted_idx = ak.argsort(events.Electron.pt, axis=1, ascending=False)
    # filter unselected electron indices
    selected_electron_idx = sorted_idx[default_mask[sorted_idx]]
    


    return events, SelectionResult(
        objects = {
            "Electron": {
                "Electron": selected_electron_idx
            }
        }
    )

@selector(
    uses=({
        f"Muon.{var}" for var in [
            # four momenta information
            "pt", "eta", "phi", "mass", 
            # quality criteria
            "isGlobal", "isStandalone", "isTracker", "nStations", "nTrackerLayers",
            # impact parameters
            "dxy", "dz", "sip3d",
            # IDs
            "tightId", "mvaId", "highPtId",
            #isolation
            "pfRelIso03_all",
        ]}
    ),
    exposed=False,
)
def muon_selector(self: Selector, events: ak.Array, **kwargs) -> tuple[ak.Array, SelectionResult]:
    min_pt = 15
    max_eta = 2.4
    selected_muon_mask = (
        # Global or Tracker Muon
        (events.Muon.isGlobal | events.Muon.isTracker) &
        # Discard Standalone Muon tracks if reconstructed in muon system only
        (~events.Muon.isStandalone | (events.Muon.nTrackerLayers >0)) &
        # WIP: Discard muons with muonBestTrackType==2 even if they are global or tracker muons
        # --> muonBestTrackType not available?
        (events.Muon.pt > min_pt) & 
        (abs(events.Muon.eta) < max_eta) &
        (events.Muon.dxy < 0.5) & (events.Muon.dz < 1.0) &
        (abs(events.Muon.sip3d) < 4) &
        # PF muon ID if pT < 200 GeV, PF muon ID or High-pT muon ID if pT > 200 GeV
        (((events.Muon.pt > 200) & events.Muon.highPtId == 2) | (events.Muon.tightId)) &
        (events.Muon.pfRelIso03_all < 0.35)
    )
    sorted_idx = ak.argsort(events.Muon.pt, axis=1, ascending=False)
    # filter unselected muon indices
    muon_idx = sorted_idx[selected_muon_mask[sorted_idx]]
    return events, SelectionResult(
        objects= {
            "Muon": {
                "Muon": muon_idx
            }
        }
    )