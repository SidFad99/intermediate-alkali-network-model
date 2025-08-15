
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

import streamlit as st


with st.expander("About the model (maths & interpretation)", expanded=True):
    st.markdown(r"""
**Conductivity**
\[
\sigma_{\text{simple}} = S\,\frac{1+\beta H}{1+\alpha H}, \qquad
H = 1 - \sum_i p_i^2, \quad p_i = s_i / \sum_j s_j.
\]
Larger \(\beta\) tends to help (near‑miss corridors); larger \(\alpha\) tends to hinder (junction blocking).

**How the box reflects your inputs**
- More/larger strands and thicker lines for larger \(p_i\).
- If 'Tie geometry to α,β' is on: thresholds and encounter bias depend on α, β so the markers respond clearly.

**References**
- Simpson index: https://en.wikipedia.org/wiki/Diversity_index#Simpson_index
- Mixed‑alkali overview: https://doi.org/10.1016/j.jnoncrysol.2010.03.013
- Jonscher UDR: https://doi.org/10.1038/267673a0
""")
with st.sidebar:
    st.header("Model inputs")
    m = st.slider("Number of ion types (m)", 1, 4, 2, 1)
    st.caption("Enter any positive 'sizes' (relative prevalence). They will be normalised.")
    sizes = []
    default_sizes = [1.0, 5.0, 1.0, 1.0]
    for i in range(m):
        sizes.append(st.number_input(f"Size s[{i+1}]", 0.0, 1000.0, default_sizes[i], 0.1, key=f"s{i}"))
    alpha = st.slider("Blocking α", 0.0, 3.0, 1.30, 0.05)
    beta  = st.slider("Mutual β",   0.0, 3.0, 0.00, 0.05)
    S_cal = st.number_input("Scale S (optional)", 0.0, 10.0, 1.0, 0.1)

    st.markdown("---")
    st.header("3‑D box controls")
    tie_geometry = st.checkbox("Tie geometry to α,β", value=True)
    base_strands = st.slider("Base strands (scaled by sizes)", 1, 6, 3, 1)
    steps = st.slider("Steps per strand", 30, 500, 200, 10)
    step_len = st.slider("Step length", 0.01, 0.08, 0.045, 0.005)
    repel = st.slider("Self‑repel radius", 0.02, 0.25, 0.10, 0.01)
    start_spread = st.slider("Start spread (0=close together)", 0.0, 0.8, 0.20, 0.05)
    encounter_bias_base = st.slider("Encounter bias (base)", 0.0, 1.0, 0.15, 0.05)
    periodic = st.checkbox("Periodic boundary (wrap)", value=True)
    d_block_base = st.slider("Blocking distance (base)", 0.005, 0.10, 0.020, 0.002)
    d_near_base  = st.slider("Near‑miss distance (base)", 0.01, 0.20, 0.06, 0.005)
    opacity = st.slider("Strand opacity", 0.1, 1.0, 0.85, 0.05)
    seed = st.number_input("Random seed", 0, 10_000, 23, 1)

sigma, H, p = sigma_simple(sizes, alpha, beta, S=S_cal)

# Effective geometry if toggled
if tie_geometry:
    # Saturating mappings (stronger, smooth response):
    # Map α in [0,3] -> multiplier in [1, ~3]
    mult_block = 1.0 + 0.7 * (alpha / (1.0 + 0.5*alpha)) * 3.0
    # Map β in [0,3] -> multiplier in [1, ~3], and widen the gap between near and block
    mult_near  = 1.0 + 0.7 * (beta  / (1.0 + 0.5*beta)) * 3.0
    d_block_eff = d_block_base * mult_block
    d_near_eff  = max(d_block_eff + 0.01, d_near_base * mult_near)
    # Encounter bias rises with β, minimal when β=0
    encounter_bias_eff = min(1.0, encounter_bias_base + 0.6 * (beta / (1.0 + beta)))
else:
    d_block_eff = d_block_base
    d_near_eff  = max(d_block_base + 0.01, d_near_base)
    encounter_bias_eff = encounter_bias_base

# Strands and linewidth scale with sizes
p_arr = np.array(p)
if p_arr.sum() == 0:
    p_scaled = np.ones(m) / m
else:
    p_scaled = p_arr / (p_arr.max() + 1e-12)
strands_each = np.maximum(1, np.round(base_strands * (0.6*p_scaled + 0.4)).astype(int))
linew_each   = 1.0 + 2.0 * p_scaled  # line width

# Generate paths per ion
paths = []
seed_base = int(seed)
for ion_idx in range(m):
    ion_paths = generate_network_paths_for_ion(
        strands_per_ion=int(strands_each[ion_idx]), steps_per=int(steps), step_len=float(step_len),
        repel=float(repel), start_spread=float(start_spread),
        encounter_bias_eff=float(encounter_bias_eff), periodic=bool(periodic),
        seed_offset=seed_base + ion_idx*97
    )
    paths.append(ion_paths)

# Analyse interactions
block_pts, mutual_pts = analyse_interactions(paths, d_block_eff=float(d_block_eff), d_near_eff=float(d_near_eff))

# ---------- Display ----------
c1, c2 = st.columns([1,1])
with c1:
    st.subheader("Simple conductivity")
    st.write(f"**σ_simple = {sigma:.4f}** (units of S)")
    st.write(f"Mixing index: **H = {H:.3f}**")
    st.write(f"Normalised sizes p: {np.array(p)}")
with c2:
    st.subheader("Visualiser mapping (live)")
    st.write(f"Effective thresholds: d_block_eff = {d_block_eff:.4f}, d_near_eff = {d_near_eff:.4f}")
    st.write(f"Encounter bias (effective): {encounter_bias_eff:.3f}")
    st.write(f"Strands per ion: {list(map(int, strands_each))}")
    st.write(f"Line widths per ion: {[round(x,2) for x in linew_each]}")

st.markdown("---")
st.subheader("3‑D simulation box (responsive to sizes and α/β)")

# Fixed palette: up to 4 consistent colours
palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

fig = plt.figure(figsize=(7.8, 6.6))
ax = fig.add_subplot(111, projection='3d')

# Box edges
for s,e in [
    ((0,0,0),(1,0,0)),((0,1,0),(1,1,0)),((0,0,1),(1,0,1)),((0,1,1),(1,1,1)),
    ((0,0,0),(0,1,0)),((1,0,0),(1,1,0)),((0,0,1),(0,1,1)),((1,0,1),(1,1,1)),
    ((0,0,0),(0,0,1)),((1,0,0),(1,0,1)),((0,1,0),(0,1,1)),((1,1,0),(1,1,1))
]:
    s=np.array(s); e=np.array(e); ax.plot([s[0],e[0]],[s[1],e[1]],[s[2],e[2]], lw=0.5, alpha=0.25, color="#999999")

# Plot paths with per-ion linewidth and a single colour per ion
for ion_idx, ion_strands in enumerate(paths):
    colour = palette[ion_idx % len(palette)]
    for s_idx, P in enumerate(ion_strands):
        ax.plot(P[:,0], P[:,1], P[:,2], lw=float(linew_each[ion_idx]), color=colour, alpha=float(opacity),
                label=(f"Ion {ion_idx+1}" if s_idx==0 else None))

# Points
if mutual_pts.shape[0]>0:
    ax.scatter(mutual_pts[:,0], mutual_pts[:,1], mutual_pts[:,2], s=8, marker='o', alpha=0.55, color="#1f78b4", label="Mutual corridors")
if block_pts.shape[0]>0:
    ax.scatter(block_pts[:,0], block_pts[:,1], block_pts[:,2], s=18, marker='^', alpha=0.9, color="#e66101", label="Blocking nodes")

ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_zlim(0,1)
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
ax.set_title(f"σ_simple={sigma:.3f}, H={H:.3f} (α={alpha}, β={beta}) — m={m}")
ax.legend(loc="upper left")
st.pyplot(fig)

st.caption("Model: σ_simple = S·(1+βH)/(1+αH),  H = 1 − ∑ p_i²,  p_i = s_i / Σ s_j. "
           "With 'Tie geometry to α,β' on, the marker density responds to α (blocking) and β (mutual).")