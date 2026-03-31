#!/usr/bin/env python3
"""
Read a LUND output file from the FARM simulation and plot
the scattered-electron kinematics and DIS variables.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# ── Constants ──────────────────────────────────────────────────────────
MP = 0.938272  # proton mass [GeV]
ME = 0.000511  # electron mass [GeV]


# ── LUND reader ────────────────────────────────────────────────────────
def read_lund(filename):
    """Parse a LUND file and return header info + particle arrays."""
    events = []
    with open(filename) as f:
        while True:
            header = f.readline()
            if not header:
                break
            hfields = header.split()
            if len(hfields) < 10:
                break
            npart = int(hfields[0])
            beam_energy = float(hfields[6])

            particles = []
            for _ in range(npart):
                pline = f.readline().split()
                particles.append({
                    "idx":    int(pline[0]),
                    "charge": int(pline[1]),
                    "pid":    int(pline[3]),
                    "mother": int(pline[4]),
                    "px":     float(pline[6]),
                    "py":     float(pline[7]),
                    "pz":     float(pline[8]),
                    "E":      float(pline[9]),
                    "mass":   float(pline[10]),
                })
            events.append({"beam_energy": beam_energy, "particles": particles})
    return events


# ── Compute DIS variables ──────────────────────────────────────────────
def compute_dis(events):
    """Extract scattered-electron kinematics and DIS variables."""
    data = dict(
        E_e=[], p_e=[], theta_e=[], phi_e=[],
        Q2=[], xB=[], nu=[], W=[], y=[]
    )

    for ev in events:
        e_beam = ev["beam_energy"]
        # Find the scattered electron (pid == 11)
        electron = None
        for p in ev["particles"]:
            if p["pid"] == 11:
                electron = p
                break
        if electron is None:
            continue

        px, py, pz, E = electron["px"], electron["py"], electron["pz"], electron["E"]
        p_mag = np.sqrt(px**2 + py**2 + pz**2)
        theta = np.arccos(pz / p_mag) if p_mag > 0 else 0.0
        phi = np.arctan2(py, px)

        # 4-momentum transfer: q = k - k'  (target at rest)
        nu = e_beam - E
        Q2 = 2.0 * e_beam * E * (1.0 - np.cos(theta))
        xB = Q2 / (2.0 * MP * nu) if nu > 0 else 0.0
        W2 = MP**2 + 2.0 * MP * nu - Q2
        W = np.sqrt(max(W2, 0.0))
        y = nu / e_beam if e_beam > 0 else 0.0

        data["E_e"].append(E)
        data["p_e"].append(p_mag)
        data["theta_e"].append(np.degrees(theta))
        data["phi_e"].append(np.degrees(phi))
        data["Q2"].append(Q2)
        data["xB"].append(xB)
        data["nu"].append(nu)
        data["W"].append(W)
        data["y"].append(y)

    return {k: np.array(v) for k, v in data.items()}


# ── Plotting ───────────────────────────────────────────────────────────
def plot(data, output_pdf="lund_plots.pdf"):
    fig, axes = plt.subplots(3, 3, figsize=(14, 12))
    fig.suptitle("Scattered Electron & DIS Variables", fontsize=16)

    hist_kw = dict(bins=100, histtype="stepfilled", alpha=0.7, color="steelblue", edgecolor="black", linewidth=0.3)

    ax = axes[0, 0]
    ax.hist(data["E_e"], **hist_kw)
    ax.set_xlabel(r"$E'_e$ [GeV]")
    ax.set_ylabel("Counts")
    ax.set_title("Electron energy")

    ax = axes[0, 1]
    ax.hist(data["p_e"], **hist_kw)
    ax.set_xlabel(r"$|\vec{p}_e|$ [GeV]")
    ax.set_title("Electron momentum")

    ax = axes[0, 2]
    ax.hist(data["theta_e"], **hist_kw)
    ax.set_xlabel(r"$\theta_e$ [deg]")
    ax.set_title("Electron polar angle")

    ax = axes[1, 0]
    ax.hist(data["phi_e"], **hist_kw)
    ax.set_xlabel(r"$\phi_e$ [deg]")
    ax.set_title("Electron azimuthal angle")

    ax = axes[1, 1]
    ax.hist(data["Q2"], **hist_kw)
    ax.set_xlabel(r"$Q^2$ [GeV$^2$]")
    ax.set_title(r"$Q^2$")

    ax = axes[1, 2]
    ax.hist(data["xB"], range=(0, 1), **hist_kw)
    ax.set_xlabel(r"$x_B$")
    ax.set_title("Bjorken x")

    ax = axes[2, 0]
    ax.hist(data["nu"], **hist_kw)
    ax.set_xlabel(r"$\nu$ [GeV]")
    ax.set_title("Energy transfer")

    ax = axes[2, 1]
    ax.hist(data["W"], **hist_kw)
    ax.set_xlabel(r"$W$ [GeV]")
    ax.set_title("Invariant mass W")

    ax = axes[2, 2]
    ax.hist(data["y"], range=(0, 1), **hist_kw)
    ax.set_xlabel("y")
    ax.set_title("Inelasticity y")

    for a in axes.flat:
        a.tick_params(labelsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_pdf, dpi=150)
    print(f"Saved → {output_pdf}")


# ── 2-D correlation plots ─────────────────────────────────────────────
def plot_2d(data, output_pdf="lund_plots_2d.pdf"):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("DIS Correlations", fontsize=16)

    hist2d_kw = dict(bins=100, cmap="inferno", cmin=1)

    ax = axes[0, 0]
    ax.hist2d(data["xB"], data["Q2"], range=[[0, 1], [0, np.percentile(data["Q2"], 99)]], **hist2d_kw)
    ax.set_xlabel(r"$x_B$"); ax.set_ylabel(r"$Q^2$ [GeV$^2$]")
    ax.set_title(r"$Q^2$ vs $x_B$")

    ax = axes[0, 1]
    ax.hist2d(data["W"], data["Q2"], range=[[0, np.percentile(data["W"], 99)], [0, np.percentile(data["Q2"], 99)]], **hist2d_kw)
    ax.set_xlabel(r"$W$ [GeV]"); ax.set_ylabel(r"$Q^2$ [GeV$^2$]")
    ax.set_title(r"$Q^2$ vs $W$")

    ax = axes[1, 0]
    ax.hist2d(data["theta_e"], data["p_e"],
              range=[[0, np.percentile(data["theta_e"], 99)], [0, np.percentile(data["p_e"], 99)]], **hist2d_kw)
    ax.set_xlabel(r"$\theta_e$ [deg]"); ax.set_ylabel(r"$|\vec{p}_e|$ [GeV]")
    ax.set_title(r"$p_e$ vs $\theta_e$")

    ax = axes[1, 1]
    ax.hist2d(data["y"], data["xB"], range=[[0, 1], [0, 1]], **hist2d_kw)
    ax.set_xlabel("y"); ax.set_ylabel(r"$x_B$")
    ax.set_title(r"$x_B$ vs $y$")

    for a in axes.flat:
        a.tick_params(labelsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_pdf, dpi=150)
    print(f"Saved → {output_pdf}")


# ── Main ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <lund_file> [output_prefix]")
        sys.exit(1)

    lund_file = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else "lund_plots"

    print(f"Reading {lund_file} ...")
    events = read_lund(lund_file)
    print(f"  {len(events)} events loaded")

    data = compute_dis(events)
    print(f"  {len(data['Q2'])} events with a scattered electron")

    plot(data, output_pdf=f"{prefix}_1d.pdf")
    plot_2d(data, output_pdf=f"{prefix}_2d.pdf")
