"""
Correction for leptons with regards to FSR photons
"""
from functools import partial

from columnflow.calibration import Calibrator, calibrator
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, Route, EMPTY_FLOAT
from columnflow.production.util import attach_coffea_behavior

from h4l.selection.lepton import electron_selection, muon_selector


np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = partial(set_ak_column, value_type=np.float32)


@calibrator(
    uses=(
        # list of external TaskArrayFunctions
        {electron_selection, muon_selector, attach_coffea_behavior} |
        {f"Photon.{var}" for var in [
            # four momenta information
            "pt", "eta", "phi", "mass", "charge",
            # isolation
            "pfRelIso03_all",
        ]} | {
            f"{field}.{var}"
            for field in ("Electron", "Muon")
            for var in [
                "pt", "eta", "phi", "mass", "charge",
            ]
        }
    ),
    produces={attach_coffea_behavior},
    exposed=False,

)
def fsr_photon_calibrator(
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
            "Muon": {
                "type_name": "Muon",
            },
            "Photon": {
                "type_name": "Photon",
            },
        },
    )
    events, ele_selection_result = self[electron_selection](
        events, working_point="loose", **kwargs,
    )
    events, muon_selection_result = self[muon_selector](
        events, working_point="loose", **kwargs,
    )
    ele_mask = ele_selection_result.x.mask
    ele = ak.mask(events.Electron, ele_mask)

    muon_mask = muon_selection_result.x.mask
    muons = ak.mask(events.Muon, muon_mask)

    # TODO: do electron/muon selections at loose working point first
    # do fsr correction for electrons
    # note: probably have to do electrons+ muons in one go to account for double
    # counting from testing in IPython:
    kinematic_photon_mask = (
        (events.Photon.pt > 2) &
        (abs(events.Photon.eta) < 2.4) &
        (events.Photon.pfRelIso03_all < 1.8)
    )
    kinematic_photons = events.Photon[kinematic_photon_mask]
    # metric_table does not exist for concatenated array, so first do the
    # DR computation per lepton flavor
    dr_ele, (ele_comb, photon_ele_comb) = ele.metric_table(
        kinematic_photons,
        return_combinations=True,
    )
    dr_mu, (mu_comb, photon_mu_comb) = muons.metric_table(
        kinematic_photons,
        return_combinations=True,
    )

    # now stitch them together at axis 1
    dr = ak.concatenate((dr_ele, dr_mu), axis=1)
    photon_comb = ak.concatenate((photon_ele_comb, photon_mu_comb), axis=1)
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

    slimmed_clean_photon_idx = ak.drop_none(clean_photon_idx)

    # One issue arises here: Photons loose their 'Photon' characteristics
    # with the sum operation. In 2/1500 cases, this leads to nan masses
    # tried re-attaching the behavior, didn't work so for. Manually
    # reset FSR calibration in these cases for now
    fsr_photons = photon_comb[slimmed_clean_photon_idx].sum(axis=-1)
    real_photon_mask = ak.num(slimmed_clean_photon_idx, axis=-1) > 0
    # finally, update the lepton four momenta
    # first, we need to figure out where the electrons and where the muons are
    # Since electrons are always first in the concatenation, we can create a mask
    # using the lepton indices
    ele_fsr_photon_mask = ak.local_index(fsr_photons) < ak.num(ele)
    updated_ele = ele.add(fsr_photons[ele_fsr_photon_mask])

    # everything else corresponds to muons
    updated_muons = muons.add(fsr_photons[~ele_fsr_photon_mask])

    # finally, correct the lepton four momenta
    for field in ["Electron", "Muon"]:
        for var in ["pt", "eta", "phi", "mass"]:
            colname = f"{field}.fsr_uncorrected_{var}"
            name = f"{field}.{var}"
            events = set_ak_column(events, colname, Route(name).apply(events, EMPTY_FLOAT))

    # The fsr photons have None entries when no photon can be matched
    # While this is physical, it leads to the side effect that the updated
    # particles (particularly electrons) can have nan masses.
    # Therefore, fill values (pt, eta, phi, mass) depending on whether there
    # actually is a photon for an electron or a muon. This is defined in the
    # real photon mask, so zip all information needed to do an ak.where
    for fieldname, new_array, backup_array, mask in zip(
        ["Electron", "Muon"],
        [updated_ele, updated_muons],
        [events.Electron, events.Muon],
        [real_photon_mask[ele_fsr_photon_mask], real_photon_mask[~ele_fsr_photon_mask]],
    ):
        masses = ak.where(mask, new_array.mass, backup_array.mass)
        nan_mask = ~np.isfinite(masses)
        if ak.any(nan_mask):
            print(f"warning: Will reset calibration for {ak.sum(nan_mask)} {fieldname}!")
            masses = ak.where(nan_mask, backup_array.mass, masses)
            # propagate the reset also to the other observables

        events = set_ak_column(events, f"{fieldname}.mass", masses)
        pts = ak.where(mask, new_array.pt, backup_array.pt)
        pts = ak.where(nan_mask, backup_array.pt, pts)
        events = set_ak_column(events, f"{fieldname}.pt", pts)

        etas = ak.where(mask, new_array.eta, backup_array.eta)
        etas = ak.where(nan_mask, backup_array.eta, etas)
        events = set_ak_column(events, f"{fieldname}.eta", etas)

        phis = ak.where(mask, new_array.phi, backup_array.phi)
        phis = ak.where(nan_mask, backup_array.phi, phis)
        events = set_ak_column(events, f"{fieldname}.phi", phis)

    return events


@fsr_photon_calibrator.init
def fsr_photon_calibrator_init(self: Calibrator) -> None:
    # add shifted new variables variables
    self.produces |= {
        f"{field}.{var}"
        for field in ("Electron", "Muon")
        for var in [
            f"{prefix}{obs}"
            for prefix in ("", "fsr_uncorrected_")
            for obs in ("pt", "eta", "phi", "mass")
        ]
    } | {"run", "luminosityBlock", "event"}
