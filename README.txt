Simple IAN Model + 3‑D Visualiser (Streamlit)
=============================================

This app shows a minimal mathematical model for total ionic conductivity with
an interactive 3‑D visualiser that highlights *mutual corridors* and *blocking nodes*.

Files
-----
- `ian_simple_gui_app.py`   — the Streamlit app (v5 engine)
- `requirements.txt`        — Python deps
- `README.txt`              — this guide

How to run
----------
1) Python 3.10+
2) Create and activate a virtual environment.
   Windows:
       python -m venv .venv
       .venv\Scripts\activate
   macOS / Linux:
       python -m venv .venv
       source .venv/bin/activate
3) Install packages:
       pip install -r requirements.txt
4) Launch the GUI:
       streamlit run ian_simple_gui_app.py


Mathematical model (simplified)
-------------------------------
You choose the number of ions `m` and positive “sizes” s₁,…,s_m (relative prevalence).
They are normalised to fractions pᵢ = sᵢ / Σⱼ sⱼ.

We define a mixing index (Simpson‑type):
    H = 1 − Σᵢ pᵢ²
H is near 0 when one ion dominates; it rises when sizes are comparable.

The total conductivity is
    σ_total = S · (1 + β H) / (1 + α H)

- α ≥ 0 is the *blocking* parameter (junctions/intersections that hinder motion).
- β ≥ 0 is the *mutualism* parameter (near‑parallel corridors that help motion).
- S is an optional scale to match units.

The form is bounded, monotone (β raises σ, α lowers σ), and smooth in the inputs.
It is intended as a compact teaching and exploration tool.


How the visualiser reflects the maths
-------------------------------------
- **Ion sizes (sᵢ):** control pᵢ and therefore H in the equation. They also scale
  the *number of strands* and the *line thickness* for each ion.
- **Blocking (α):** increases the blocking detection radius d_block,eff → more orange triangles.
  In the equation, α increases the denominator → σ falls.
- **Mutualism (β):** increases the near‑miss detection radius d_near,eff and a small
  encounter bias so strands come closer → more blue dots. In the equation, β increases
  the numerator → σ rises.
- A sidebar toggle (“Tie geometry to α,β”) lets you switch this coupling on/off for demos.

Tips (especially for m = 2)
---------------------------
- Use multiple strands per ion, small start spread, periodic boundaries, and a modest
  encounter bias to ensure the strands meet.
- Keep `blocking distance` < `near‑miss distance`. Raise both to see more markers.

References (trusted)
--------------------
- Simpson diversity index (definition & properties):
  https://en.wikipedia.org/wiki/Diversity_index#Simpson_index
- Review context on mixed‑alkali effects and pathway interactions:
  https://doi.org/10.1016/j.jnoncrysol.2010.03.013
- Universal trends for ionic conduction in disordered media (background):
  Jonscher, *Nature* 267, 673 (1977) — https://doi.org/10.1038/267673a0
