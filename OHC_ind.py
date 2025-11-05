"""
Python port of MATLAB OHC_ind.m (Verhulst 2018 model helper)
- Flexible audiogram frequency lists (any order)
- Robust path handling (spaces ok)
- Outputs Poles/<name>/profile.txt and StartingPoles.dat
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# Try both loaders for .mat files
try:
    from scipy.io import loadmat
except ImportError:
    loadmat = None

try:
    import h5py
except ImportError:
    h5py = None

# ----------------------------- Utils -----------------------------

def load_mat_var(path, candidates):
    """Load a .mat file and return the first matching variable by name."""
    if isinstance(candidates, str):
        candidates = [candidates]
    
    # Try h5py first (for v7.3 files)
    if h5py is not None:
        try:
            with h5py.File(path, 'r') as f:
                keys = list(f.keys())
                for c in candidates:
                    if c in f:
                        data = f[c][:]
                        # If it's a reference to another dataset, dereference it
                        while isinstance(data, np.ndarray) and data.size == 1:
                            ref = data.flat[0]
                            if isinstance(ref, np.ndarray) and ref.dtype == 'object':
                                data = f[ref.flat[0]][:]
                            else:
                                break
                        return data
                if len(keys) == 1:
                    data = f[keys[0]][:]
                    while isinstance(data, np.ndarray) and data.size == 1:
                        ref = data.flat[0]
                        if isinstance(ref, np.ndarray) and ref.dtype == 'object':
                            data = f[ref.flat[0]][:]
                        else:
                            break
                    return data
        except (OSError, ValueError):
            # Not an HDF5 file, fall through to scipy
            pass
    
    # Fallback to scipy's loadmat
    if loadmat is not None:
        M = loadmat(path)
        # Remove meta-keys
        keys = [k for k in M.keys() if not k.startswith('__')]
        for c in candidates:
            if c in M:
                return M[c]
        # If only one non-meta key, return it
        if len(keys) == 1:
            return M[keys[0]]
    
    raise KeyError(f"None of the expected variables {candidates} found in {os.path.basename(path)}")

def ensure_1d(arr):
    arr = np.asarray(arr)
    return np.ravel(arr)

def nearest_index(array, value):
    array = np.asarray(array)
    return int(np.argmin(np.abs(array - value)))

def write_lines(path, lines):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        for line in lines:
            f.write(line.rstrip() + '\n')

def semilogx_like(ax, x, y, **kwargs):
    ax.plot(x, y, **kwargs)
    ax.set_xscale('log')

# ------------------------ Core computation -----------------------

def ohc_ind(
    name,
    hl_freqs_hz,
    hl_db,
    base_dir=".",
    show_figs=True,
):
    """
    Parameters
    ----------
    name : str
        Output profile folder name under Poles/<name>.
    hl_freqs_hz : list[float]
        Frequencies (Hz) for the provided audiogram points (any order).
    hl_db : list[float]
        Hearing loss values (dB) corresponding to hl_freqs_hz.
    base_dir : str
        Directory containing 'mat files', 'Poles/Flat00/StartingPoles.dat', etc.
    show_figs : bool
        Whether to display matplotlib figures.
    """

    # --- Paths (handle spaces safely) ---
    mat_dir = os.path.join(base_dir, "mat files")
    poles_flat_dir = os.path.join(base_dir, "Poles", "Flat00")
    out_dir = os.path.join(base_dir, "Poles", name)

    # --- Load model data ---
    # Present for completeness; not directly used here, but loaded to mirror MATLAB
    _ = load_mat_var(os.path.join(mat_dir, "BWrange.mat"), ["BWrange", "BWrange_Hz"])  # unused downstream

    cf = ensure_1d(load_mat_var(os.path.join(mat_dir, "cf.mat"), "cf")).astype(float)  # Hz
    ModelQ = ensure_1d(load_mat_var(os.path.join(mat_dir, "ModelQ.mat"), ["ModelQ", "Q"])).astype(float)

    PT = loadmat(os.path.join(mat_dir, "PoleTrajs.mat"))
    # Try to fetch by common names; adapt if your file uses different names
    BMS = PT.get("BMS", None)     # cell/array of spectra over Poles
    SMax = PT.get("SMax", None)   # typically [nPoles x 1000]
    SN   = PT.get("SN", None)     # often abscissa (pole index vs CF index)
    if SMax is None:
        raise KeyError("Could not find SMax in PoleTrajs.mat. Please check variable names.")
    SMax = np.asarray(SMax)
    if SN is not None:
        SN = np.asarray(SN)

    # Powerlaw parameters for Q computation
    PL = loadmat(os.path.join(mat_dir, "Powerlawpar.mat"))
    # Look for 'a' and 'b' (scalars or 1x1 arrays)
    a = float(np.ravel(PL.get("a", np.array([np.nan])))[0])
    b = float(np.ravel(PL.get("b", np.array([np.nan])))[0])
    if not np.isfinite(a) or not np.isfinite(b):
        raise KeyError("Powerlawpar.mat must contain scalars 'a' and 'b'.")

    # StartingPoles (Flat00)
    sp_path = os.path.join(poles_flat_dir, "StartingPoles.dat")
    SP_all = np.loadtxt(sp_path, dtype=float)
    # remove first entry (middle ear)
    StartingPoles = ensure_1d(SP_all)[1:]  # length should be 1000
    NHSP = StartingPoles.copy()

    # Poles axis
    Poles = np.round(np.arange(0.036, 0.302 + 1e-12, 0.001), 3)  # 0.036:0.001:0.302

    # --- Find nearest NH pole index for each CF section (choose closest above/below) ---
    NHn = np.empty_like(NHSP, dtype=int)
    for k in range(len(NHSP)):
        # first index where Poles > NHSP[k]
        greater = np.where(Poles > NHSP[k])[0]
        if greater.size == 0:
            # All Poles <= NHSP[k]; choose last
            NHn[k] = len(Poles) - 1
            continue
        Hn = greater[0]
        if Hn == 0:
            NHn[k] = 0
        else:
            Ln = Hn - 1
            # pick closer
            if (Poles[Hn] - NHSP[k]) < (NHSP[k] - Poles[Ln]):
                NHn[k] = Hn
            else:
                NHn[k] = Ln

    # --- GainDiff over poles relative to NH starting pole for each CF ---
    n_poles = len(Poles)
    n_cf = len(NHSP)
    if SMax.shape[0] != n_poles or SMax.shape[1] != n_cf:
        raise ValueError(
            f"SMax shape {SMax.shape} does not match expected (nPoles={n_poles}, nCF={n_cf})."
        )
    GainDiff = (SMax[NHn, np.arange(n_cf)] - SMax)  # broadcast over rows

    # ---------------- Audiogram → HL over CF sections ----------------
    hl_freqs_hz = np.asarray(hl_freqs_hz, dtype=float)
    hl_db = np.asarray(hl_db, dtype=float)
    if hl_freqs_hz.shape != hl_db.shape:
        raise ValueError("hl_freqs_hz and hl_db must have the same length.")

    # Map each provided freq to nearest CF index
    probe_idx = [nearest_index(cf, f) for f in hl_freqs_hz]
    # Deduplicate indices by keeping the last provided HL for any duplicate index
    probe_map = {}
    for idx, hl in zip(probe_idx, hl_db):
        probe_map[idx] = float(hl)
    # Sort by ascending index along the cochlea
    probes_sorted = np.array(sorted(probe_map.items()), dtype=object)  # [[idx, hl], ...]
    probes = probes_sorted[:, 0].astype(int)
    hls_at_probes = probes_sorted[:, 1].astype(float)

    # Build HL array over all CF sections via piecewise linear interpolation on index axis
    HL_sections = np.zeros_like(cf, dtype=float)

    # fill before first probe
    HL_sections[:probes[0]] = hls_at_probes[0]
    # linear segments
    for i in range(len(probes) - 1):
        i0, i1 = probes[i], probes[i + 1]
        y0, y1 = hls_at_probes[i], hls_at_probes[i + 1]
        if i1 > i0:
            x = np.arange(i0, i1 + 1, dtype=int)
            # linear interpolation on index axis
            HL_sections[i0:i1 + 1] = y0 + (y1 - y0) * (x - i0) / max(1, (i1 - i0))
        else:
            # Should not happen after sorting; safeguard
            HL_sections[i1:i0 + 1] = np.linspace(y1, y0, i0 - i1 + 1)
    # fill after last probe
    HL_sections[probes[-1]:] = hls_at_probes[-1]

    # ---------------- Compute HI poles via threshold crossing ----------------
    HISPn = np.empty(n_cf, dtype=int)
    for k in range(n_cf):
        # first n where GainDiff[n,k] > HL_sections[k]
        n_found = np.where(GainDiff[:, k] > HL_sections[k])[0]
        if n_found.size == 0:
            HISPn[k] = n_poles - 1  # fallback (last valid index, equivalent to MATLAB's 267)
        else:
            HISPn[k] = int(n_found[0])
    HISP = Poles[HISPn]

    # Prepend one pole for middle ear section (to match MATLAB convention)
    StartingPolesHI = np.concatenate([[HISP[0]], HISP])
    StartingPolesNH = np.concatenate([[NHSP[0]], NHSP])

    # ModelQHI from power law
    ModelQHI = a * np.power(StartingPolesHI, b)

    # ----------------------------- Plots -----------------------------
    # 1) Pole trajectories (if BMS/SN present)
    if BMS is not None and SN is not None:
        fig1, ax1 = plt.subplots()
        # Emulate the MATLAB loop (plot every 100th)
        for n in range(0, 1000, 100):
            # BMS may be cell-like; try to index safely
            try:
                # If BMS is an object array/cell-like
                bn = BMS.ravel()[n]
                y = np.ravel(bn)
                ax1.plot(y)
                ax1.plot(np.ravel(SN[:, n]), np.ravel(SMax[:, n]), 'rs', markersize=3)
            except Exception:
                # Fallback: skip plotting BMS if format unknown
                pass
        ax1.set_xlim(0, 2000)
        ax1.set_title('The Pole Trajectories for Each CF')
        if show_figs:
            plt.show()

    # 2) Audiogram shaping
    fig2, ax2 = plt.subplots()
    semilogx_like(ax2, cf, HL_sections, color='r', linewidth=2, label='Interpolated HL')
    # Scatter the provided audiogram points
    ax2.semilogx(hl_freqs_hz, hl_db, 'bo', linewidth=2, label='Audiogram points')
    ax2.set_xlim(125, 25000)
    ax2.set_ylim(-10, 50)
    ax2.set_yticks(np.arange(-10, 55, 10))
    ax2.invert_yaxis()
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Hearing Loss [dB]')
    ax2.legend(loc='best', frameon=False)
    if show_figs:
        plt.show()

    # 3) Starting poles & Q
    fig3, (ax3a, ax3b) = plt.subplots(2, 1, figsize=(5, 10))
    ax3a.plot(StartingPoles, 'k-', linewidth=1, label='NHfit')
    ax3a.plot(StartingPolesHI[1:], 'r-', linewidth=2, label='HIfit')  # drop ME to match NH length
    ax3a.set_ylim(0, 0.3)
    ax3a.set_xlabel('cochlear section')
    ax3a.set_ylabel('Starting Pole Value')
    ax3a.legend(loc='upper right', frameon=False)

    # Q vs CF (prepend cf[0] like MATLAB)
    cf_with_me = np.concatenate([[cf[0]], cf])
    ax3b.semilogx(cf_with_me / 1000.0, ModelQ, 'b-', label='NH')
    ax3b.semilogx(cf_with_me / 1000.0, ModelQHI, 'r-', label='HI')
    ax3b.set_xlim(0.04, 23)
    ax3b.set_ylim(0, 20)
    ax3b.set_xlabel('CF [kHz]')
    ax3b.set_ylabel(r'$Q_{ERB}$')
    ax3b.legend(loc='upper left', frameon=False)
    plt.tight_layout()
    if show_figs:
        plt.show()

    # ----------------------------- Save outputs -----------------------------
    os.makedirs(out_dir, exist_ok=True)

    # profile.txt: line1=cf, line2=HL (section-wise)
    prof_path = os.path.join(out_dir, "profile.txt")
    line1 = "\t".join(f"{x:.6f}" for x in cf)
    line2 = "\t".join(f"{x:.6f}" for x in HL_sections)
    write_lines(prof_path, [line1, line2])

    # StartingPoles.dat
    sp_out_path = os.path.join(out_dir, "StartingPoles.dat")
    if name == "Normal":
        np.savetxt(sp_out_path, StartingPolesNH, fmt="%.6E")
    else:
        np.savetxt(sp_out_path, StartingPolesHI, fmt="%.6E")

    return {
        "cf": cf,
        "HL_sections": HL_sections,
        "StartingPolesNH": StartingPolesNH,
        "StartingPolesHI": StartingPolesHI,
        "ModelQ": ModelQ,
        "ModelQHI": ModelQHI,
        "output_dir": out_dir,
    }

# Example usage (commented out to prevent execution on import)
# ohc_ind(name='Brent', 
#         hl_freqs_hz=[8000, 6000, 4000, 3000, 2000, 1000, 500, 250, 125], 
#         hl_db=[0, 0, 0, 0, 0, 0, 0, 0, 0])