#!/usr/bin/env python3
"""
Compare the Hertz3D simulation (sphere on flat, third-medium contact) with the
analytical Hertz solution.

Parameters are NOT hardcoded: they are read from the files that define the
simulation, so they cannot drift apart:
  * hertz3D_test.geo -> R, void_gap, Hf, W, arc_ratio
  * Hertz3D.cxx      -> H (applied displacement), te, E, nu

Usage:
  python hertz_validation.py                       # defaults below
  python hertz_validation.py --force out/force.txt --geo hertz3D_test.geo \
                             --cxx Hertz3D.cxx --component 1
Any parameter can be overridden on the command line (--R, --gap, --E, ...).

Model facts used (quarter model, symmetry planes x=0 and z=0):
  * the force measured on the model is 1/4 of the physical force -> x4
  * the imposed displacement u(t) = H t/te starts with a gap, so contact only
    starts at u = gap and the Hertz approach is  delta = u - gap
  * identical elastic bodies, sphere on a plane:  R* = R,
    1/E* = 2 (1 - nu^2) / E
"""
import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np


# --------------------------------------------------------------------------
# Reading the parameters from the .geo / .cxx files
# --------------------------------------------------------------------------
def _strip_cpp_comments(text):
    return re.sub(r"//.*", "", text)


def parse_geo(path):
    """Evaluate the simple `name = expression;` assignments of a gmsh file."""
    env = {"Pi": math.pi, "Sin": math.sin, "Cos": math.cos,
           "Max": max, "Min": min, "Sqrt": math.sqrt}
    for line in _strip_cpp_comments(Path(path).read_text()).splitlines():
        m = re.match(r"\s*(\w+)\s*=\s*([^;{}]+);", line)
        if not m:
            continue
        try:
            env[m.group(1)] = eval(m.group(2), {"__builtins__": {}}, env)
        except Exception:
            pass  # loop-dependent or non-numeric lines: not needed
    return env


def parse_cxx(path):
    """Read `constexpr auto X = mfem_mgis::real{value};` constants."""
    text = _strip_cpp_comments(Path(path).read_text())
    pat = (r"constexpr\s+(?:const\s+)?auto\s+(\w+)\s*=\s*"
           r"mfem_mgis::real\s*\{\s*([^}]+?)\s*\}")
    return {m.group(1): float(m.group(2)) for m in re.finditer(pat, text)}


# --------------------------------------------------------------------------
# Hertz solution (sphere / plane)
# --------------------------------------------------------------------------
def effective_modulus(E1, nu1, E2, nu2):
    return 1.0 / ((1 - nu1**2) / E1 + (1 - nu2**2) / E2)


def hertz_force(delta, R, Estar):
    delta = np.clip(delta, 0.0, None)
    return 4.0 / 3.0 * Estar * math.sqrt(R) * delta**1.5


def hertz_contact_radius(delta, R):
    return np.sqrt(R * np.clip(delta, 0.0, None))


def hertz_delta_from_force(F, R, Estar):
    return (3.0 * np.clip(F, 0.0, None) / (4.0 * Estar * math.sqrt(R))) ** (2.0 / 3.0)


# --------------------------------------------------------------------------
# Reading the force file written by ComputeResultantForceOnBoundary
# --------------------------------------------------------------------------
def load_force(path, component=None):
    rows = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            rows.append([float(v) for v in line.replace(",", " ").split()])
        except ValueError:
            continue  # header line
    if not rows:
        sys.exit(f"No numeric data found in {path}")
    data = np.array(rows)
    t = data[:, 0]
    forces = data[:, 1:]
    if forces.shape[1] == 1:
        return t, forces[:, 0], 0
    if component is None:
        # the contact force is the dominant component (loading is along y)
        component = int(np.argmax(np.abs(forces[-1])))
        print(f"[info] force file has {forces.shape[1]} components, "
              f"using component {component} (largest). Override: --component")
    return t, forces[:, component], component


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", default="output/force.txt", help="force file from the simulation")
    ap.add_argument("--geo", default="hertz3D_test.geo")
    ap.add_argument("--cxx", default="Hertz3D.cxx")
    ap.add_argument("--component", type=int, default=None,
                    help="column of the force vector (0=x,1=y,2=z); auto if omitted")
    ap.add_argument("--symmetry-factor", type=float, default=4.0,
                    help="physical force / force in the model (quarter model -> 4)")
    ap.add_argument("--out", default="hertz_validation.png")
    # optional overrides
    for k in ("R", "gap", "E", "nu", "H", "te"):
        ap.add_argument(f"--{k}", type=float, default=None)
    ap.add_argument("--no-show", action="store_true")
    a = ap.parse_args()

    geo = parse_geo(a.geo) if Path(a.geo).exists() else {}
    cxx = parse_cxx(a.cxx) if Path(a.cxx).exists() else {}

    def pick(cli, src, key, name):
        if cli is not None:
            return cli
        if key in src:
            return src[key]
        sys.exit(f"Cannot determine '{name}': pass --{name if name != 'gap' else 'gap'}")

    R = pick(a.R, geo, "R", "R")
    gap = pick(a.gap, geo, "void_gap", "gap")
    E = pick(a.E, cxx, "E", "E")
    nu = pick(a.nu, cxx, "nu", "nu")
    H = pick(a.H, cxx, "H", "H")
    te = pick(a.te, cxx, "te", "te")
    r_tm = geo.get("R", R) * math.sin(geo.get("arc_ratio", 0.5) * math.pi / 2) \
        if "arc_ratio" in geo else None

    Estar = effective_modulus(E, nu, E, nu)
    print(f"R={R}  gap={gap}  E={E:.3e}  nu={nu}  H={H}  te={te}  E*={Estar:.4e}")

    t, f_model, comp = load_force(a.force, a.component)
    F_fem = np.abs(f_model) * a.symmetry_factor
    u = H * t / te                      # imposed displacement of IND_TOP
    delta = u - gap                     # Hertz approach (contact starts at u = gap)

    F_h = hertz_force(delta, R, Estar)
    a_h = hertz_contact_radius(delta, R)
    # indentation the Hertz law would need to carry the FEM force: the
    # difference with `delta` is the bulk compliance of the finite bodies
    delta_eq = hertz_delta_from_force(F_fem, R, Estar)

    # ---- sanity checks on the validity of the analytic reference ----------
    d_max = delta.max()
    print(f"max indentation delta/R = {d_max / R:.3f}   "
          f"max contact radius a/R = {a_h.max() / R:.3f}")
    if d_max / R > 0.1:
        print("[warn] delta/R > 0.1: finite-strain regime, Hertz (small strain, "
              "half-space) is only a first-order reference.")
    if r_tm is not None and a_h.max() > r_tm:
        print(f"[warn] Hertz contact radius {a_h.max():.3f} exceeds the third-medium "
              f"patch radius {r_tm:.3f}: contact leaves the meshed void region.")

    # ---- table -------------------------------------------------------------
    mask = delta > 0.02 * max(d_max, 1e-30)
    print(f"\n{'t':>8} {'u':>9} {'delta':>9} {'F_fem':>12} {'F_hertz':>12} {'ratio':>7}")
    idx = np.linspace(0, len(t) - 1, min(10, len(t))).astype(int)
    for i in idx:
        r = F_fem[i] / F_h[i] if F_h[i] > 0 else float("nan")
        print(f"{t[i]:8.3f} {u[i]:9.4f} {delta[i]:9.4f} {F_fem[i]:12.4e} {F_h[i]:12.4e} {r:7.3f}")
    if mask.any():
        err = np.abs(F_fem[mask] / F_h[mask] - 1)
        print(f"\nmean |F_fem/F_hertz - 1| over contact steps: {100 * err.mean():.1f} %")

    # ---- plots -------------------------------------------------------------
    import matplotlib
    if a.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))

    dd = np.linspace(0, max(d_max, 1e-12), 300)
    ax[0].plot(dd / R, hertz_force(dd, R, Estar), "k-", label="Hertz")
    ax[0].plot(np.clip(delta, 0, None) / R, F_fem, "o", ms=4, label="FEM (third medium)")
    ax[0].axvline(0, color="grey", lw=0.5)
    ax[0].set_xlabel(r"$\delta / R = (u - \mathrm{gap})/R$")
    ax[0].set_ylabel("Total contact force F")
    ax[0].set_title("Force vs indentation")
    ax[0].legend()
    ax[0].grid(alpha=0.3)

    ax[1].plot(u / R, F_fem, "o", ms=4, label="FEM")
    ax[1].plot(u / R, F_h, "k-", label="Hertz")
    ax[1].axvline(gap / R, color="r", ls="--", lw=1, label="gap closes")
    ax[1].set_xlabel("Imposed displacement u / R")
    ax[1].set_ylabel("F")
    ax[1].set_title("Force vs imposed displacement")
    ax[1].legend()
    ax[1].grid(alpha=0.3)

    ax[2].plot(delta[mask] / R, F_fem[mask] / F_h[mask], "o-", ms=4)
    ax[2].axhline(1, color="k", lw=1)
    ax[2].set_xlabel(r"$\delta / R$")
    ax[2].set_ylabel(r"$F_{FEM} / F_{Hertz}$")
    ax[2].set_title("Ratio (1 = perfect agreement)")
    ax[2].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print(f"\nFigure written to {a.out}")
    if not a.no_show:
        plt.show()


if __name__ == "__main__":
    main()