"""
Correction for leptons with regards to FSR photons
"""
from functools import partial

from columnflow.calibration import Calibrator, calibrator
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, Route, EMPTY_FLOAT
from columnflow.production.util import attach_coffea_behavior

from h4l.selection.lepton import electron_selection


np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = partial(set_ak_column, value_type=np.float32)


@calibrator(
    uses=(
        # list of external TaskArrayFunctions
        {electron_selection, attach_coffea_behavior} |
        {f"Photon.{var}" for var in [
            # four momenta information
            "pt", "eta", "phi", "mass", "charge",
            # isolation
            "pfRelIso03_all",
        ]} | {
            f"Electron.{var}" for var in [
                "pt", "eta", "phi", "mass", "charge",
            ]
        }
    ),
    produces=set(),
    exposed=False,

)
def electron_fsr_photon_calibrator(
    self: Calibrator,
    events: ak.Array,
    check_double_asignment: bool = True,
    **kwargs,
):
    # first, attach coffea behavior for four momenta logic
    events = self[attach_coffea_behavior](
        events,
        collections={
            "Electron": {
                "type_name": "Electron",
            },
            "Photon": {
                "type_name": "Photon",
            },
        },
    )
    events, selection_result = self[electron_selection](
        events, working_point="loose", **kwargs,
    )
    ele_mask = selection_result.x.mask
    ele = ak.mask(events.Electron, ele_mask)
    # TODO: do electron/muon selections at loose working point first
    # do fsr correction for electrons
    # note: probably have to do electrons+ muons in one go to account for double
    # counting from testing in IPython:
    kinematic_photon_mask = (
        (events.Photon.pt > 2) &
        (abs(events.Photon.eta) < 2.4) &
        (events.Photon.pfRelIso03_all < 1.8)
    )
    dr, (ele_comb, photon_comb) = ele.metric_table(
        events.Photon[kinematic_photon_mask],
        return_combinations=True,
    )
    # get indices of photon combinations to check for double counting

    # apply critieria, e.g.
    weighed_dr = dr / photon_comb.energy**2
    mask = (dr < 0.5) & (weighed_dr < 0.12)
    photon_idx = ak.local_index(photon_comb)
    if check_double_asignment:
        def find_minimum_idx(
            current_indices=None,
            mask=None,
            max_iterations=4,
            current_iterations=0,
        ):
            if current_iterations >= max_iterations:
                raise ValueError(
                    f"Was unable to match all photons to leptons after {max_iterations} iterations!" # noqa
                )
            if mask is None:
                mask = ak.ones_like(weighed_dr, type=bool)
            reduced_weighed_dr = ak.mask(weighed_dr, mask)
            # first, identify the minimum value in each column
            col_minval = ak.min(reduced_weighed_dr, axis=1, keepdims=True)

            # now broadcast this to the same shape as the original array
            # note that ak.broadcast_arrays returns a list with the broadcasted
            # array at index 0 and the original array at index 1
            broadcast_col_minval = ak.broadcast_arrays(col_minval, reduced_weighed_dr)[0]

            # the minimum positions are now the exactly where the original array
            # matches the broadcasted array with minimum values
            min_mask = reduced_weighed_dr == broadcast_col_minval

            # the minimum indices are the photon indices where the above mask is True
            min_idx = photon_idx[min_mask]
            # add the photon indices with the minimal weighed dr to the list
            # of matched photon indices
            # from IPython import embed; embed()

            if current_indices is None:
                final_indices = min_idx
            else:
                final_indices = ak.concatenate((current_indices, min_idx), axis=-1)

            # since we only match photons with the minimum weighted dr
            # it's still possible to have unmatched photons, so look for them

            # if a photon has already been matched to a lepton, there is an entry
            # in the corresponding column. Use ak.any to check for this
            matched_idx_mask = ak.any(min_mask, axis=1, keepdims=True)

            # identify unmatched photons by reversing the selection
            unmatched_idx_mask = ~matched_idx_mask

            # if there are still any unmatched photons that would in principle
            # fullfil our selection criteria, continue
            newmask = (mask & unmatched_idx_mask)
            if ak.any(newmask):
                return find_minimum_idx(
                    current_indices=final_indices,
                    mask=newmask,
                    max_iterations=max_iterations,
                    current_iterations=current_iterations + 1,
                )
            return final_indices

        clean_photon_idx = find_minimum_idx(mask=mask, max_iterations=ak.max(photon_idx))

    else:
        clean_photon_idx = photon_idx[mask]

    fsr_photons = photon_comb[ak.drop_none(clean_photon_idx)].sum(axis=-1)

    # finally, update the lepton four momenta
    updated_ele = ele.add(fsr_photons)

    # finally, correct the lepton four momenta
    for var in ["pt", "eta", "phi", "mass"]:
        colname = f"Electron.fsr_uncorrected_{var}"
        name = f"Electron.{var}"
        events = set_ak_column(events, colname, Route(name).apply(events, EMPTY_FLOAT))

    for fieldname, array in zip(["Electron"], [updated_ele]):
        events = set_ak_column(events, f"{fieldname}.pt", array.pt)
        events = set_ak_column(events, f"{fieldname}.eta", array.eta)
        events = set_ak_column(events, f"{fieldname}.phi", array.phi)
        events = set_ak_column(events, f"{fieldname}.mass", array.mass)

    return events
