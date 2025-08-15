
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# --- Self-contained simplified model (no external import needed) ---
def sigma_simple(sizes, alpha, beta, S=1.0):
    sizes = [float(x) for x in sizes]
    tot = sum(sizes)
    if tot <= 0:
        return 0.0, 0.0, [0.0]*len(sizes)
    p = [x/tot for x in sizes]
    H = 1.0 - sum(pi*pi for pi in p)
    sigma = S * (1.0 + beta*H) / (1.0 + alpha*H)
    return sigma, H, p

st.set_page_config(page_title="Simple IAN Model + 3‑D Visualiser", layout="wide")

with st.expander("About the model (maths & interpretation)", expanded=True):
    st.write("Conductivity: sigma_simple = S * (1 + beta * H) / (1 + alpha * H).")
    st.write("H = 1 - sum(p_i^2), with p_i = s_i / sum_j s_j.")
    st.write("beta increases conductivity (mutual corridors); alpha decreases it (blocking).")
    st.write("The box reflects your inputs: more/larger strands for larger p_i; thresholds and encounter bias respond to alpha and beta.")
    st.write("References: Simpson index (Wikipedia), Mixed-alkali overview (J. Non-Cryst. Solids 2010), Jonscher UDR (Nature 1977).")

# ---------- Geometry helpers ----------
def _rand_unit_vec(rng):
    v = rng.normal(size=3); n = np.linalg.norm(v) + 1e-12
    return v / n

def _wrap_periodic(p): return p - np.floor(p)

def generate_network_paths_for_ion(strands_per_ion, steps_per, step_len, repel, 
                                   start_spread, encounter_bias_eff, periodic, seed_offset):
    rng = np.random.default_rng(seed_offset)
    ion_strands = []; centre = np.array([0.5, 0.5, 0.5])
    for _ in range(strands_per_ion):
        start = centre + (rng.uniform(-1,1,size=3) * start_spread * 0.5)
        start = np.clip(start, 0, 1); pts = [start]; v = _rand_unit_vec(rng)
        for _ in range(steps_per):
            p = pts[-1]
            bounce = np.zeros(3)
            if not periodic:
                for d in range(3):
                    if p[d] < 0.05: bounce[d] += (0.08 - p[d])
                    if p[d] > 0.95: bounce[d] -= (p[d] - 0.92)
            repel_vec = np.zeros(3)
            for q in pts[-12:]:
                r = p - q; r2 = np.dot(r,r)
                if r2 > 1e-6 and r2 < repel**2:
                    repel_vec += r / (r2 + 1e-6)
            pull = (centre - p)
            dirn = 0.70*v + 0.55*bounce + 0.20*repel_vec + encounter_bias_eff*0.25*pull
            nrm = np.linalg.norm(dirn)
            dirn = _rand_unit_vec(rng) if nrm < 1e-10 else (dirn / nrm)
            newp = p + step_len*dirn; newp = _wrap_periodic(newp) if periodic else np.clip(newp, 0, 1)
            pts.append(newp); v = 0.6*v + 0.4*dirn; v = v / (np.linalg.norm(v) + 1e-12)
        ion_strands.append(np.vstack(pts))
    return ion_strands

def _segment_distance(a0, a1, b0, b1):
    A = a1 - a0; B = b1 - b0; C = a0 - b0
    AA = np.dot(A, A) + 1e-12; BB = np.dot(B, B) + 1e-12
    AB = np.dot(A, B); AC = np.dot(A, C); BC = np.dot(B, C)
    den = AA*BB - AB*AB + 1e-18
    s = np.clip((AB*BC - BB*AC)/den, 0.0, 1.0); t = np.clip((AA*BC - AB*AC)/den, 0.0, 1.0)
    p = a0 + s*A; q = b0 + t*B
    return np.linalg.norm(p - q), p, q

def analyse_interactions(paths, d_block_eff, d_near_eff, max_points=(2000, 3000)):
    block_pts = []; mutual_pts = []; m = len(paths)
    for i in range(m):
        for j in range(i+1, m):
            for Pi in paths[i]:
                for Pj in paths[j]:
                    for k in range(len(Pi)-1):
                        for l in range(len(Pj)-1):
                            d, p, q = _segment_distance(Pi[k], Pi[k+1], Pj[l], Pj[l+1])
                            if d < d_block_eff: block_pts.append(0.5*(p+q))
                            elif d < d_near_eff: mutual_pts.append(0.5*(p+q))
    import numpy as np
    if len(block_pts) > max_points[0]:
        idx = np.random.choice(len(block_pts), max_points[0], replace=False)
        block_pts = [block_pts[k] for k in idx]
    if len(mutual_pts) > max_points[1]:
        idx = np.random.choice(len(mutual_pts), max_points[1], replace=False)
        mutual_pts = [mutual_pts[k] for k in idx]
    return np.array(block_pts) if block_pts else np.zeros((0,3)), np.array(mutual_pts) if mutual_pts else np.zeros((0,3))

# ---------- Sidebar ----------
st.sidebar.header("Model inputs")
m = st.sidebar.slider("Number of ion types (m)", 1, 4, 2, 1)
st.sidebar.caption("Enter any positive 'sizes' (relative prevalence). They will be normalised.")
sizes = [st.sidebar.number_input(f"Size s[{i+1}]", 0.0, 1000.0, 1.0 if i == 0 else 5.0 if i == 1 else 1.0, 0.1, key=f"s{i}") for i in range(m)]
alpha = st.sidebar.slider("Blocking α", 0.0, 3.0, 1.30, 0.05)
beta  = st.sidebar.slider("Mutual β",   0.0, 3.0, 0.00, 0.05)
S_cal = st.sidebar.number_input("Scale S (optional)", 0.0, 10.0, 1.0, 0.1)

st.sidebar.markdown("---")
st.sidebar.header("3‑D box controls")
tie_geometry = st.sidebar.checkbox("Tie geometry to α,β", value=True)
base_strands = st.sidebar.slider("Base strands (scaled by sizes)", 1, 6, 3, 1)
steps = st.sidebar.slider("Steps per strand", 30, 500, 200, 10)
step_len = st.sidebar.slider("Step length", 0.01, 0.08, 0.045, 0.005)
repel = st.sidebar.slider("Self‑repel radius", 0.02, 0.25, 0.10, 0.01)
start_spread = st.sidebar.slider("Start spread (0=close together)", 0.0, 0.8, 0.20, 0.05)
encounter_bias_base = st.sidebar.slider("Encounter bias (base)", 0.0, 1.0, 0.15, 0.05)
periodic = st.sidebar.checkbox("Periodic boundary (wrap)", value=True)
d_block_base = st.sidebar.slider("Blocking distance (base)", 0.005, 0.10, 0.020, 0.002)
d_near_base  = st.sidebar.slider("Near‑miss distance (base)", 0.01, 0.20, 0.06, 0.005)
opacity = st.sidebar.slider("Strand opacity", 0.1, 1.0, 0.85, 0.05)
seed = st.sidebar.number_input("Random seed", 0, 10_000, 23, 1)

# ---------- Model ----------
sigma, H, p = sigma_simple(sizes, alpha, beta, S=S_cal)

if tie_geometry:
    mult_block = 1.0 + 0.7 * (alpha / (1.0 + 0.5*alpha)) * 3.0
    mult_near  = 1.0 + 0.7 * (beta  / (1.0 + 0.5*beta)) * 3.0
    d_block_eff = d_block_base * mult_block
    d_near_eff  = max(d_block_eff + 0.01, d_near_base * mult_near)
    encounter_bias_eff = min(1.0, encounter_bias_base + 0.6 * (beta / (1.0 + beta)))
else:
    d_block_eff = d_block_base
    d_near_eff  = max(d_block_base + 0.01, d_near_base)
    encounter_bias_eff = encounter_bias_base

p_arr = np.array(p); p_scaled = (p_arr / (p_arr.max() + 1e-12)) if p_arr.sum() > 0 else (np.ones(m)/m)
strands_each = np.maximum(1, np.round(base_strands * (0.6*p_scaled + 0.4)).astype(int))
linew_each   = 1.0 + 2.0 * p_scaled

# Generate & analyse
paths = []
seed_base = int(seed)
for ion_idx in range(m):
    ion_paths = generate_network_paths_for_ion(
        strands_per_ion=int(strands_each[ion_idx]), steps_per=int(steps), step_len=float(step_len),
        repel=float(repel), start_spread=float(start_spread),
        encounter_bias_eff=float(encounter_bias_eff), periodic=bool(periodic),
        seed_offset=seed_base + ion_idx*97
    ); paths.append(ion_paths)

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
palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

fig = plt.figure(figsize=(7.8, 6.6)); ax = fig.add_subplot(111, projection='3d')
# Box edges
for s,e in [((0,0,0),(1,0,0)),((0,1,0),(1,1,0)),((0,0,1),(1,0,1)),((0,1,1),(1,1,1)),
            ((0,0,0),(0,1,0)),((1,0,0),(1,1,0)),((0,0,1),(0,1,1)),((1,0,1),(1,1,1)),
            ((0,0,0),(0,0,1)),((1,0,0),(1,0,1)),((0,1,0),(0,1,1)),((1,1,0),(1,1,1))]:
    s=np.array(s); e=np.array(e); ax.plot([s[0],e[0]],[s[1],e[1]],[s[2],e[2]], lw=0.5, alpha=0.25, color="#999999")

for ion_idx, ion_strands in enumerate(paths):
    colour = palette[ion_idx % len(palette)]
    for s_idx, P in enumerate(ion_strands):
        ax.plot(P[:,0], P[:,1], P[:,2], lw=float(linew_each[ion_idx]), color=colour, alpha=float(opacity),
                label=(f"Ion {ion_idx+1}" if s_idx==0 else None))

if mutual_pts.shape[0]>0:
    ax.scatter(mutual_pts[:,0], mutual_pts[:,1], mutual_pts[:,2], s=8, marker='o', alpha=0.55, color="#1f78b4", label="Mutual corridors")
if block_pts.shape[0]>0:
    ax.scatter(block_pts[:,0], block_pts[:,1], block_pts[:,2], s=18, marker='^', alpha=0.9, color="#e66101", label="Blocking nodes")

ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_zlim(0,1)
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
ax.set_title(f"σ_simple={sigma:.3f}, H={H:.3f} (α={alpha}, β={beta}) — m={m}")
ax.legend(loc="upper left"); st.pyplot(fig)

st.caption("Model: sigma_simple = S*(1+beta*H)/(1+alpha*H),  H = 1 − sum(p_i^2),  p_i = s_i / sum(s_j).")
