#!/usr/bin/env python3
"""
Shared-error reanalysis - main figures.
Fig. 1  GPT-4 vs human crowd (primary archive):   a scatter, b same-way share, c direction split
Fig. 3  GPT-4 vs domain experts (secondary archive): same template
Fig. 2  Each predictor's error vs GPT-4's (primary archive), small multiples

These are the figures of the Matters Arising manuscript (fig1-3.png in
output/figures/). The same code is embedded, eval-guarded, in
docs/analysis_report.Rmd; R/05_figures.R draws simpler R diagnostic
versions of Figs 1-2 (F1_*.png, F2_*.png), not the manuscript figures.

Run from anywhere:  python3 scripts/make_figures.py
Inputs are read from the repository's data/derived/ and outputs land in
output/figures/, both resolved relative to this script's location.
Requires: pandas, numpy, matplotlib
"""
import os, math
import pandas as pd, numpy as np, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
from statistics import NormalDist
trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))

# ============================ CONFIG ============================
REPO           = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ALL_MODELS_CSV = os.path.join(REPO, "data", "derived", "comparisons_all_models.csv")
ARCHIVE2_CSV   = os.path.join(REPO, "data", "derived", "comparisons_archive2.csv")
OUTDIR         = os.environ.get("FIG_OUTDIR", os.path.join(REPO, "output", "figures"))
BOOT           = 2000        # bootstrap reps for observed-rate / slope CIs
BOOT_NULL      = 400         # bootstrap reps for null CIs (each recomputes the integral)
# ================================================================

os.makedirs(OUTDIR, exist_ok=True)
N = NormalDist()
T = 1.96; CLIP = 12; VIEW = 13.4

mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":9.5,"axes.edgecolor":"#9AA1A9",
    "axes.linewidth":0.8,"xtick.color":"#4A525B","ytick.color":"#4A525B",
    "text.color":"#222222","axes.labelcolor":"#222222"})

RED   = "#A93B2B"   # same-side points / observed bars
BLUE  = "#4C6FA1"   # opposite-side points
LBLUE = "#A9BFD8"   # opposite-direction bar
GRAYP = "#C9CCD1"   # non-joint-miss points
GRAY1 = "#ABABAB"   # corrected-null bar
GRAY2 = "#D7D7D7"   # naive-null bar
INK   = "#222222"
FILL  = "#F6ECEA"   # quadrant shading

clip = lambda v: np.clip(v, -CLIP, CLIP)
_erf = np.vectorize(math.erf)
def Phi(x): return 0.5*(1.0+_erf(np.asarray(x,dtype=float)/1.4142135623730951))

def null_exact(miss_g, miss_h, t=T, ngrid=1201):
    """Corrected null (vectorized): joint-miss rate, same-side rate, direction share
    for two independent predictors scored against the same noisy estimate."""
    if miss_g <= 0 or miss_h <= 0: return {"joint":0.0,"same":0.0,"direction":0.5}
    sig = lambda p: t/N.inv_cdf(1-p/2)
    s2g = max(sig(miss_g)**2-1, 0); s2h = max(sig(miss_h)**2-1, 0)
    u = np.linspace(-9, 9, ngrid); phi = np.exp(-u*u/2)/math.sqrt(2*math.pi)
    def hi(s2): return 1-Phi((t+u)/math.sqrt(s2)) if s2 > 0 else (-u > t).astype(float)
    def lo(s2): return Phi((-t+u)/math.sqrt(s2)) if s2 > 0 else (-u < -t).astype(float)
    gh, gl, hh, hl = hi(s2g), lo(s2g), hi(s2h), lo(s2h)
    joint = trapz((gh+gl)*(hh+hl)*phi, u)
    same  = trapz((gh*hh+gl*hl)*phi, u)
    return {"joint":joint, "same":same, "direction":same/joint}

def cluster_groups(c):
    uq = np.unique(c); return [np.where(c == k)[0] for k in uq]

def slope_ci(x, y, c, B=BOOT):
    np.random.seed(42)
    b = ((x-x.mean())*(y-y.mean())).sum()/((x-x.mean())**2).sum()
    groups = cluster_groups(c); K = len(groups); est = []
    for _ in range(B):
        idx = np.concatenate([groups[j] for j in np.random.randint(0, K, K)])
        xi, yi = x[idx], y[idx]; d = ((xi-xi.mean())**2).sum()
        if d > 0: est.append(((xi-xi.mean())*(yi-yi.mean())).sum()/d)
    return b, *np.percentile(est, [2.5, 97.5])

def boot_stats(x, y, c, B=BOOT, Bn=BOOT_NULL):
    """Cluster bootstrap CIs for: same-way share (of all), direction share (of joint),
    and the two null benchmarks in panel b."""
    np.random.seed(42)
    groups = cluster_groups(c); K = len(groups)
    sw, dr, cn, nn = [], [], [], []
    for i in range(B):
        idx = np.concatenate([groups[j] for j in np.random.randint(0, K, K)])
        xg, yg = x[idx], y[idx]
        gm = np.abs(xg) > T; hm = np.abs(yg) > T; both = gm & hm
        same = both & (np.sign(xg) == np.sign(yg))
        sw.append(same.mean())
        if both.sum() > 0: dr.append(same.sum()/both.sum())
        if i < Bn and gm.mean() > 0 and hm.mean() > 0:
            ne = null_exact(gm.mean(), hm.mean(), ngrid=601)
            cn.append(ne["same"]); nn.append(gm.mean()*hm.mean()*0.5)
    pct = lambda a: (100*np.percentile(a, 2.5), 100*np.percentile(a, 97.5))
    return pct(sw), pct(dr), pct(cn), pct(nn)

# ---------------- the Fig-1 template (used for Fig. 1 and Fig. 3) ----------------
def panel_fig(x, y, c, other_name, unit_word, title_a, outpath):
    n = len(x)
    gm = np.abs(x) > T; hm = np.abs(y) > T; both = gm & hm
    same = both & (np.sign(x) == np.sign(y)); oppo = both & ~same
    over  = same & (x > T); under = same & (x < -T)
    n_joint, n_same, n_oppo = int(both.sum()), int(same.sum()), int(oppo.sum())
    ne = null_exact(gm.mean(), hm.mean())
    obs_b  = 100*same.mean();  corr_b = 100*ne["same"];  naive_b = 100*gm.mean()*hm.mean()*0.5
    obs_c  = 100*n_same/n_joint; oppo_c = 100*n_oppo/n_joint; corr_c = 100*ne["direction"]
    ci_sw, ci_dr, ci_cn, ci_nn = boot_stats(x, y, c)
    ci_op = (100-ci_dr[1], 100-ci_dr[0])

    fig = plt.figure(figsize=(11.7, 6.6), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.5, 1], height_ratios=[1, 1])
    axa = fig.add_subplot(gs[:, 0]); axb = fig.add_subplot(gs[0, 1]); axc = fig.add_subplot(gs[1, 1])

    # ---- a: scatter ----
    axa.add_patch(Rectangle((T, T), VIEW-T, VIEW-T, facecolor=FILL, edgecolor="none", zorder=0))
    axa.add_patch(Rectangle((-VIEW, -VIEW), VIEW-T, VIEW-T, facecolor=FILL, edgecolor="none", zorder=0))
    axa.axhline(0, color="#B9BDC2", lw=.8, zorder=1); axa.axvline(0, color="#B9BDC2", lw=.8, zorder=1)
    for s in (T, -T):
        axa.axhline(s, color="#7A828B", lw=.8, ls=(0, (4, 3)), zorder=1)
        axa.axvline(s, color="#7A828B", lw=.8, ls=(0, (4, 3)), zorder=1)
    oth = ~both
    axa.scatter(clip(x[oth]),  clip(y[oth]),  s=5,  c=GRAYP, alpha=.6, lw=0, zorder=2)
    axa.scatter(clip(x[oppo]), clip(y[oppo]), s=11, c=BLUE,  alpha=.9, lw=0, zorder=3)
    axa.scatter(clip(x[same]), clip(y[same]), s=11, c=RED,   alpha=.9, lw=0, zorder=4)
    b, lo, hi = slope_ci(x, y, c); a0 = y.mean()-b*x.mean()
    xx = np.linspace(-CLIP, CLIP, 50)
    axa.fill_between(xx, a0+lo*xx, a0+hi*xx, color="#666666", alpha=.18, lw=0, zorder=5)
    axa.plot(xx, a0+b*xx, color="#1A1A1A", lw=1.6, zorder=6)
    axa.text(-VIEW+0.7,  VIEW-0.8, f"opposite directions\n{n_oppo}",  color=BLUE, fontsize=9, va="top")
    axa.text( VIEW-0.7,  VIEW-0.8, f"both overshoot\n{int(over.sum())}",  color=RED, fontsize=9, va="top", ha="right")
    axa.text(-VIEW+0.7, -VIEW+0.8, f"both undershoot\n{int(under.sum())}", color=RED, fontsize=9, va="bottom")
    axa.text( VIEW-0.7, -VIEW+0.8, f"slope {b:.2f}\n[{lo:.2f}, {hi:.2f}]", color=INK, fontsize=9, va="bottom", ha="right")
    axa.set_xlim(-VIEW, VIEW); axa.set_ylim(-VIEW, VIEW); axa.set_aspect("equal")
    axa.set_xticks([-10, -5, 0, 5, 10]); axa.set_yticks([-10, -5, 0, 5, 10])
    axa.set_xlabel("GPT-4 forecasting error\n(standard errors of the observed effect)", fontsize=9.5)
    axa.set_ylabel(f"{other_name} forecasting error\n(standard errors of the observed effect)", fontsize=9.5)
    axa.set_title(title_a, fontsize=10.5, pad=10)
    for sp in ("top", "right"): axa.spines[sp].set_visible(False)
    axa.text(-0.16, 1.03, "a", transform=axa.transAxes, fontsize=14, fontweight="bold")

    # ---- b: same-way share of all ----
    xs = [0, 1, 2]; vals = [obs_b, corr_b, naive_b]; cols = [RED, GRAY1, GRAY2]
    cis = [ (obs_b, *ci_sw), (corr_b, *ci_cn), (naive_b, *ci_nn) ]
    axb.grid(axis="y", color="#DDDDDD", lw=.7, zorder=0)
    axb.bar(xs, vals, width=.6, color=cols, zorder=2)
    for xpos, (v, l, h) in zip(xs, cis):
        axb.errorbar(xpos, v, yerr=[[max(v-l, 0)], [max(h-v, 0)]], fmt="none",
                     ecolor=INK, elinewidth=1.0, capsize=3, zorder=3)
        axb.text(xpos, h+0.9, f"{v:.1f}%", ha="center", fontsize=9.5, fontweight="bold")
    axb.set_xticks(xs)
    axb.set_xticklabels(["Observed", "Corrected null", "Null ignoring\nshared estimate"], fontsize=8.5)
    top = max(ci_sw[1], obs_b)*1.35
    axb.set_ylim(0, top); axb.set_yticks([0, 10, 20] if top < 30 else [0, 10, 20, 30])
    axb.set_yticklabels([f"{t:.0f}%" if t else "0" for t in axb.get_yticks()])
    axb.set_title(f"Share of all {n:,} {unit_word}\nthat both predictors miss the same way",
                  fontsize=10, loc="left", pad=8)
    for sp in ("top", "right"): axb.spines[sp].set_visible(False)
    axb.text(-0.14, 1.12, "b", transform=axb.transAxes, fontsize=14, fontweight="bold")

    # ---- c: direction split among joint misses ----
    axc.grid(axis="y", color="#DDDDDD", lw=.7, zorder=0)
    axc.bar([0, 1], [obs_c, oppo_c], width=.55, color=[RED, LBLUE], zorder=2)
    for xpos, v, (l, h) in [(0, obs_c, ci_dr), (1, oppo_c, ci_op)]:
        axc.errorbar(xpos, v, yerr=[[max(v-l, 0)], [max(h-v, 0)]], fmt="none",
                     ecolor=INK, elinewidth=1.0, capsize=3, zorder=3)
        axc.text(xpos, h+3.5, f"{v:.0f}%", ha="center", fontsize=9.5, fontweight="bold")
    axc.axhline(corr_c, color=INK, lw=1.2, ls=(0, (5, 2)), zorder=4)
    axc.text(1.42, corr_c+2, f"corrected null, {corr_c:.0f}%", fontsize=8.5, ha="right", color=INK)
    axc.axhline(50, color="#9AA1A9", lw=1.0, ls=(0, (1, 2)), zorder=4)
    axc.text(1.42, 52, "chance, 50%", fontsize=8.5, ha="right", color="#9AA1A9")
    axc.set_xticks([0, 1]); axc.set_xticklabels(["Same\ndirection", "Opposite\ndirections"], fontsize=8.5)
    axc.set_xlim(-.55, 1.55); axc.set_ylim(0, 118); axc.set_yticks([0, 50, 100])
    axc.set_yticklabels(["0", "50%", "100%"])
    axc.set_title(f"Of the {n_joint} {unit_word} both miss,\nthe share missed in the same direction",
                  fontsize=10, loc="left", pad=8)
    for sp in ("top", "right"): axc.spines[sp].set_visible(False)
    axc.text(-0.14, 1.12, "c", transform=axc.transAxes, fontsize=14, fontweight="bold")

    fig.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"{os.path.basename(outpath)}: slope {b:.2f} [{lo:.2f},{hi:.2f}] | same-way {obs_b:.1f}% "
          f"(nulls {corr_b:.1f}/{naive_b:.1f}) | direction {obs_c:.0f}% (null {corr_c:.0f})")

# ---------------- Fig 2: small-multiple scatters (unchanged style) ----------------
def scatter_small(ax, x, y, c, title, xlab):
    gm = np.abs(x) > T; em = np.abs(y) > T; both = gm & em
    same = both & (np.sign(x) == np.sign(y)); oppo = both & ~same; oth = ~both
    over = same & (x > T); under = same & (x < -T)
    ax.axhline(0, color="#DfE3E7", lw=.6); ax.axvline(0, color="#DfE3E7", lw=.6)
    for s in (T, -T):
        ax.axhline(s, color="#7A828B", lw=.6, ls=(0, (4, 3))); ax.axvline(s, color="#7A828B", lw=.6, ls=(0, (4, 3)))
    ax.scatter(clip(x[oth]),  clip(y[oth]),  s=5,  c=GRAYP, alpha=.55, lw=0)
    ax.scatter(clip(x[oppo]), clip(y[oppo]), s=7,  c=BLUE,  alpha=.85, lw=0)
    ax.scatter(clip(x[same]), clip(y[same]), s=7,  c=RED,   alpha=.85, lw=0)
    b, lo, hi = slope_ci(x, y, c, B=800); a0 = y.mean()-b*x.mean(); xx = np.linspace(-CLIP, CLIP, 40)
    ax.fill_between(xx, a0+lo*xx, a0+hi*xx, color="#666666", alpha=.16, lw=0)
    ax.plot(xx, a0+b*xx, color="#1A1A1A", lw=1.3)
    ax.text(-VIEW+0.6,  VIEW-0.6, f"opposite\n{int(oppo.sum())}", color=BLUE, fontsize=5.8, va="top", ha="left", linespacing=.9)
    ax.text( VIEW-0.6,  VIEW-0.6, f"both over\n{int(over.sum())}", color=RED, fontsize=5.8, va="top", ha="right", linespacing=.9)
    ax.text(-VIEW+0.6, -VIEW+0.6, f"both under\n{int(under.sum())}", color=RED, fontsize=5.8, va="bottom", ha="left", linespacing=.9)
    ax.text( VIEW-0.6, -VIEW+0.6, f"slope {b:.2f}\n[{lo:.2f}, {hi:.2f}]", color=INK, fontsize=5.8, va="bottom", ha="right", linespacing=.9)
    ax.set_xlim(-VIEW, VIEW); ax.set_ylim(-VIEW, VIEW); ax.set_aspect("equal")
    ax.set_xticks([-10, 0, 10]); ax.set_yticks([-10, 0, 10])
    ax.set_title(title, loc="left", fontweight="bold", fontsize=9.5)
    ax.set_ylabel("forecasting error (SEs)", fontsize=8)
    if xlab: ax.set_xlabel("GPT-4 forecasting error (SEs)", fontsize=8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)

# ============================ DATA ============================
am = pd.read_csv(ALL_MODELS_CSV)
key = ["study", "outcome.name", "reference_condition", "condition.name"]
wide = am.pivot_table(index=key, columns="model", values="z_recal", aggfunc="first").reset_index()
g = wide["gpt-4"].values; cl = wide["study"].values

# ---- Fig 1 ----
h = wide["human"].values; ok = ~(np.isnan(g) | np.isnan(h))
panel_fig(g[ok], h[ok], cl[ok], "Human", "comparisons",
          "The two predictors' errors line up, comparison by comparison",
          os.path.join(OUTDIR, "fig1.png"))

# ---- Fig 2 ----
models = [("deepseek/deepseek-chat-v3-0324", "DeepSeek-V3"), ("openai/gpt-oss-120b", "GPT-OSS-120B"),
          ("gpt-3.5-turbo", "GPT-3.5"), ("google/gemma-3-27b-it", "Gemma-3-27B"),
          ("babbage-002", "babbage-002"), ("davinci-002", "davinci-002")]
letters = "abcdef"
fig, axs = plt.subplots(2, 4, figsize=(12.4, 6.9), constrained_layout=True); axs = axs.ravel()
for i, (mc, nm) in enumerate(models):
    zo = wide[mc].values; ok = ~(np.isnan(g) | np.isnan(zo))
    scatter_small(axs[i], g[ok], zo[ok], cl[ok], f"{letters[i]}   {nm}", xlab=(i >= 2))
axs[6].axis("off"); axs[7].axis("off")
LEG = [Line2D([], [], marker='o', ls='', mfc=RED,   mec='none', label='both miss, same side'),
       Line2D([], [], marker='o', ls='', mfc=BLUE,  mec='none', label='both miss, opposite side'),
       Line2D([], [], marker='o', ls='', mfc=GRAYP, mec='none', label='not a joint miss'),
       Line2D([], [], color="#1A1A1A", lw=1.6, label='fitted slope (95% CI)')]
axs[6].legend(handles=LEG, loc="center", frameon=False, fontsize=9,
              title="each point = one of 1,678 comparisons", title_fontsize=9)
fig.savefig(os.path.join(OUTDIR, "fig2.png"), dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig2.png: done")

# ---- Fig 3 ----
a2 = pd.read_csv(ARCHIVE2_CSV).dropna(subset=["z_recal_gpt", "z_recal_expert"])
panel_fig(a2["z_recal_gpt"].values, a2["z_recal_expert"].values, a2["family"].values,
          "Expert", "conditions",
          "The two predictors' errors line up, condition by condition",
          os.path.join(OUTDIR, "fig3.png"))
print("done: fig1.png, fig2.png, fig3.png")
