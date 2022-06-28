from typing import Callable, List, Union
from cytoflow import *

import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import scipy.stats as st
from tqdm import tqdm
import scipy.optimize
import sys
from lmfit import *
from matplotlib.lines import Line2D

mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["font.family"] = "sans-serif"
plt.rcParams["figure.figsize"] = (10, 6)  # type: ignore
mpl.rcParams["pdf.fonttype"] = 42
sns.set_style(
    "ticks",
    {
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "font_color": "k",
        "axes.edgecolor": "k",
        "xtick.color": "k",
        "ytick.color": "k",
    },
)
sns.set_context("talk", font_scale=1.0)  # type: ignore


# zero-inflated-poisson PDF
def zip_pdf(
    log_cit: Union[List, float, np.ndarray],
    bprime: float,
    lmbda: float,
    sigma_on: float,
    sigma_off: float,
    mu_off: float,
) -> np.ndarray:
    # bprime = beta / (beta + alpha)
    # lmbda = k1 / gamma
    if isinstance(log_cit, list):
        lcs = log_cit
    elif isinstance(log_cit, np.ndarray):
        lcs = list(log_cit)
    else:
        lcs = [log_cit]

    p_ifoff = (1 - bprime) * st.norm(mu_off, sigma_off).pdf(log_cit)

    n_vals = np.arange(0, 1500)  # 0, 1, 2, ..., 1500
    pcgns = np.array(
        [
            np.sum(
                st.norm((n_vals / 600) + 6.5, sigma_on).pdf(c)
                * st.poisson(lmbda).pmf(n_vals)
            )
            for c in lcs
        ]
    )

    p_ifon = bprime * pcgns

    return p_ifon + p_ifoff


print("Getting data")

dx1 = pd.read_csv(
    "../02_activator_validations/data/20220303/all_cells_live_mch_gated.csv"
)
dx1["date"] = "2022-03-03"
dx2 = pd.read_csv(
    "../02_activator_validations/data/20211115/all_cells_live_mch_gated.csv"
)
dx2["date"] = "2021-11-15"
dx = pd.concat([dx1, dx2])

all_plasmids = sorted(list(set(dx["plasmid"])))
all_dates = sorted(list(set(dx["date"])))

#  [p for p in list(itertools.product(all_plasmids, all_dates)) if is_valid_pd_tuple(p*)]


def is_valid_pd_tuple(plasmid: float, date: str) -> bool:
    return dx[(dx["plasmid"] == plasmid) & (dx["date"] == date)].shape[0] > 0


pd_tuples = [
    t
    for t in list(itertools.product(all_plasmids, all_dates))
    if is_valid_pd_tuple(t[0], t[1])
]
print(len(pd_tuples))
#  sys.exit(0)
#  print(list(dx))
#  print(len(all_plasmids))
#  print(all_dates)
#  sys.exit(0)

if not os.path.exists("./popt_df.csv"):
    print("Fitting telegraph model and plotting fits")

    #  fig, ax = plt.subplots(6, 4, figsize=(16, 15))
    #  fig, ax = plt.subplots(7, 7, figsize=(19, 26))
    fig, ax = plt.subplots(10, 6, figsize=(22.5, 22.5))

    def make_plasmid_plot(plasmid: float, date: str, ax: plt.Axes) -> pd.DataFrame:
        #  print("\tMaking plot for plasmid", plasmid, "on date", date)
        df = dx[dx["plasmid"] == plasmid]
        df = df[df["date"] == date]
        df = df[df["day"] == 2]
        df = df[(df["P1"]) & (df["mCherry"])]
        if df.shape[0] > 5000:
            df = df.sample(n=4999)
        if df.shape[0] < 50:
            print(plasmid, date, df.shape)

        kde = st.gaussian_kde(np.log10(df["mCitrine-A"]))  # fit KDE to log citrine

        bp_guess = 1 - kde.integrate_box_1d(6, 7)

        name = list(df["description"])[0]

        #  print("\tFitting parameters")
        popt, _ = scipy.optimize.curve_fit(
            zip_pdf,
            np.linspace(6, 9),
            kde.pdf(np.linspace(6, 9)),
            p0=[bp_guess, 750, 0.25, 0.15, 6.3],
            bounds=[[1e-5, 1e-5, 0.1, 0.1, 6], [1, 2000, 5, 5, 7.5]],
        )
        #     print("Optimal parameters:", popt)

        log_c = np.linspace(5, 10, 500)
        y = zip_pdf(log_c, *popt)
        z = np.log10(df["mCitrine-A"])

        a2 = ax.twinx()
        sns.lineplot(x=log_c, y=y, ax=ax, color="blue")
        sns.kdeplot(x=z, ax=a2, color="red")

        for a in [ax, a2]:
            a.set_xlim(5, 10)
            a.set_xlabel("Log Citrine")
        ax.set_ylabel("Likelihood")
        a2.set_ylabel("")
        a2.set_yticks([])
        ax.set_title(name + "\n" + date)

        bprime, lmbda, sigma_on, sigma_off, mu_off = popt

        #  print("\tDone!")
        #     print(f_on, beta_guess, "\t", beta, alpha, k1, gamma, sigma)
        return pd.DataFrame.from_dict(
            {
                "plasmid": [plasmid],
                "date": [date],
                "description": [name],
                "bprime": [bprime],
                "lambda": [lmbda],
                "sigma_on": [sigma_on],
                "sigma_off": [sigma_off],
                "mu_off": [mu_off],
            }
        )

    popt_df = pd.concat(
        [
            make_plasmid_plot(pd_tuples[i][0], pd_tuples[i][1], ax.flat[i])
            for i in tqdm(range(len(pd_tuples)))
            #  pd.concat(
            #      [
            #          make_plasmid_plot(all_plasmids[i], d, ax.flat[i])
            #          # for i in tqdm(range(3))
            #          for i in tqdm(range(len(all_plasmids)))
            #      ]
            #  )
            #  for d in all_dates
        ]
    )

    custom_lines = [
        Line2D([0], [0], color="red", lw=4),
        Line2D([0], [0], color="blue", lw=4),
    ]
    ax.flat[0].legend(custom_lines, ["Real", "Fit"], loc="center right")

    sns.despine(fig)
    plt.tight_layout()

    plt.savefig("./plots/fits.pdf", bbox_inches="tight")

    popt_df.to_csv("./popt_df.csv", index=False)
else:
    print("Reading parameter dataframe from disk")
    popt_df = pd.read_csv("./popt_df.csv")

print("Plotting citrine vs beta")

val_df = pd.read_csv("../02_activator_validations/data/20220303/activation_data.csv")
#  print(list(val_df))
#  val_df.head()

pdf = (
    popt_df.set_index("plasmid")
    .join(val_df.set_index("plasmid"), on="plasmid", how="left", rsuffix="_val")
    .reset_index()
)
#  print(list(pdf))
pdf = pdf[
    [
        "plasmid",
        "description",
        "avg_enrichment_d2",
        "composition",
        "character",
        "Citrine On",
        "mCitrine-A",
        "bprime",
        "lambda",
        "sigma_on",
        "sigma_off",
        "mu_off",
    ]
]
pdf = pdf.dropna()

mfi_df = pd.read_csv("../02_activator_validations/data/mfi_data.csv")
pdf = (
    pdf.set_index("plasmid")
    .join(mfi_df.set_index("plasmid"), on="plasmid", how="left", rsuffix="_mfi")
    .reset_index()
)

fig, ax = plt.subplots(figsize=(3, 3))

sns.scatterplot(
    data=pdf,
    x="Citrine On",
    y="bprime",
    color="white",
    edgecolor="tab:blue",
    s=40,
    linewidth=2,
    ax=ax,
)

r, p = st.pearsonr(pdf["Citrine On"], pdf["bprime"])

ax.set_ylim(-0.05, 1.05)
ax.set_xlim(-0.05, 1.05)

ax.set_xticks([0, 0.5, 1])
ax.set_yticks([0, 0.5, 1])

ax.text(x=0, y=0.92, s="Pearson R$^2$={:.2f}".format(r), fontsize=16)

ax.set_xlabel("Fraction On")
ax.set_ylabel(r"$\beta'$", rotation="horizontal", ha="right")

fig.savefig("./plots/fon_vs_bprime.pdf", bbox_inches="tight")


print("Plotting citrine vs lambda")

fig, ax = plt.subplots(figsize=(3, 3))

sns.scatterplot(
    data=pdf,
    x="Citrine On",
    y="lambda",
    color="white",
    edgecolor="tab:blue",
    s=40,
    linewidth=2,
    ax=ax,
)

r, p = st.pearsonr(pdf["Citrine On"], pdf["lambda"])

ax.set_ylim(-50, 1500)
ax.set_xlim(-0.03, 1.03)

ax.set_xticks([0, 0.5, 1])
ax.set_yticks([0, 750, 1500])

ax.text(x=0, y=1325, s="Pearson R$^2$={:.2f}".format(r), fontsize=16)

ax.set_xlabel("Fraction On")
ax.set_ylabel(r"$\lambda$", rotation="horizontal", ha="right")

fig.savefig("./plots/fon_vs_lambda.pdf", bbox_inches="tight")

print("Plotting citrine vs mu_off")

fig, ax = plt.subplots(figsize=(3, 3))

sns.scatterplot(
    data=pdf,
    x="Citrine On",
    y="mu_off",
    color="white",
    edgecolor="tab:blue",
    s=40,
    linewidth=2,
    ax=ax,
)

r, p = st.pearsonr(pdf["Citrine On"], pdf["sigma_off"])

ax.set_ylim(6, 7.6)
ax.set_xlim(-0.03, 1.03)

ax.set_xticks([0, 0.5, 1])
ax.set_yticks([6, 6.5, 7, 7.5])

ax.text(x=0, y=7.43, s="Pearson R$^2$={:.2f}".format(r), fontsize=16)

ax.set_xlabel("Fraction On")
ax.set_ylabel(r"$\mu_{\mathrm{off}}$", rotation="horizontal", ha="right")

fig.savefig("./plots/fon_vs_muoff.pdf", bbox_inches="tight")

print("Plotting MFI vs lambda")

fig, ax = plt.subplots(figsize=(3, 3))

sns.scatterplot(
    data=pdf,
    x="mfi",
    y="lambda",
    color="white",
    edgecolor="tab:blue",
    s=40,
    linewidth=2,
    ax=ax,
)

r, p = st.pearsonr(pdf["mfi"], pdf["lambda"])

ax.set_ylim(-50, 1500)
ax.set_yticks([0, 500, 1000, 1500])

ax.set_xscale("log")
ax.set_xlim(1e6, 1e9)
ax.set_xticks([1e6, 1e7, 1e8, 1e9])

ax.text(x=1.05e6, y=1325, s="Pearson R$^2$={:.2f}".format(r), fontsize=16)

ax.set_xlabel("On Peak MFI")
ax.set_ylabel(r"$\lambda$", rotation="horizontal", ha="right")

fig.savefig("./plots/mfi_vs_lambda.pdf", bbox_inches="tight")

f = sns.pairplot(pdf)

f.savefig("./plots/correlogram.pdf", bbox_inches="tight")


def relu(
    x: Union[np.ndarray, List, float], intercept: float, slope: float, threshold: float
) -> float:
    return intercept + (x >= threshold) * slope * (x - threshold)  # type: ignore


print("Plotting screen score vs beta")

fig, ax = plt.subplots(figsize=(3, 3))

bp_popt, pcov = scipy.optimize.curve_fit(
    relu,
    pdf["avg_enrichment_d2"],
    pdf["bprime"],
    p0=[0, 1, 0],
    bounds=[[0, 0, -6], [1, 10, 6]],
)
xv = np.linspace(-3, 6)
yv = relu(xv, *bp_popt)
ax.plot(xv, yv, color="red", lw=2, linestyle="--", zorder=-10)

sns.scatterplot(
    data=pdf,
    x="avg_enrichment_d2",
    y="bprime",
    color="white",
    edgecolor="tab:blue",
    s=40,
    linewidth=2,
    ax=ax,
)

r, p = st.pearsonr(pdf["avg_enrichment_d2"], pdf["bprime"])

ax.set_ylim(-0.05, 1.05)
ax.set_xlim(-3, 6)

ax.set_xticks([-3, 0, 3, 6])
ax.set_yticks([0, 0.5, 1])

ax.text(x=-2.75, y=0.92, s="Pearson R$^2$={:.2f}".format(r), fontsize=16)

ax.set_xlabel("Act. $\log_2$(ON:OFF)")
ax.set_ylabel(r"$\beta'$", rotation="horizontal", ha="right")

fig.savefig("./plots/screend2_vs_bprime.pdf", bbox_inches="tight")

print("Plotting screen score vs lambda")

fig, ax = plt.subplots(figsize=(3, 3))

lm_popt, pcov = scipy.optimize.curve_fit(
    relu,
    pdf["avg_enrichment_d2"],
    pdf["lambda"],
    p0=[0, 1, 0],
    bounds=[[0, 0, -6], [1500, 15000, 6]],
)
xv = np.linspace(-3, 6, num=1000)
yv = relu(xv, *lm_popt)
ax.plot(xv, yv, color="red", lw=2, linestyle="--", zorder=-10)

sns.scatterplot(
    data=pdf,
    x="avg_enrichment_d2",
    y="lambda",
    color="white",
    edgecolor="tab:blue",
    s=40,
    linewidth=2,
    ax=ax,
)

r, p = st.pearsonr(pdf["avg_enrichment_d2"], pdf["lambda"])

ax.set_xlim(-3, 6)
ax.set_xticks([-3, 0, 3, 6])

ax.set_ylim(-50, 1500)
ax.set_yticks([0, 750, 1500])

ax.text(x=-2.75, y=1325, s="Pearson R$^2$={:.2f}".format(r), fontsize=16)

ax.set_xlabel("Act. $\log_2$(ON:OFF)")
ax.set_ylabel(r"$\lambda$", rotation="horizontal", ha="right")

fig.savefig("./plots/screend2_vs_lambda.pdf", bbox_inches="tight")

print("Plotting screen score vs mu_off")

fig, ax = plt.subplots(figsize=(3, 3))

mu_popt, pcov = scipy.optimize.curve_fit(
    relu,
    pdf["avg_enrichment_d2"],
    pdf["mu_off"],
    p0=[6.5, 1, 0],
    bounds=[[6, 0, -6], [8, 20, 6]],
)
print(mu_popt)
xv = np.linspace(-3, 6, num=1000)
yv = relu(xv, *mu_popt)
ax.plot(xv, yv, color="red", lw=2, linestyle="--", zorder=-10)

sns.scatterplot(
    data=pdf,
    x="avg_enrichment_d2",
    y="mu_off",
    color="white",
    edgecolor="tab:blue",
    s=40,
    linewidth=2,
    ax=ax,
)

r, p = st.pearsonr(pdf["avg_enrichment_d2"], pdf["mu_off"])

ax.set_ylim(6, 7.6)
ax.set_xlim(-3.05, 6.05)

ax.set_xticks([-3, 0, 3, 6])
ax.set_yticks([6, 6.5, 7, 7.5])

ax.text(x=-2.75, y=7.43, s="Pearson R$^2$={:.2f}".format(r), fontsize=16)

ax.set_xlabel("Act. $\log_2$(ON:OFF)")
ax.set_ylabel(r"$\mu_{\mathrm{off}}$", rotation="horizontal", ha="right")

fig.savefig("./plots/screen_d2_vs_muoff.pdf", bbox_inches="tight")

print("Plotting screen score vs sigma_off")

fig, ax = plt.subplots(figsize=(3, 3))

so_popt, pcov = scipy.optimize.curve_fit(
    relu,
    pdf["avg_enrichment_d2"],
    pdf["sigma_off"],
    p0=[0, 1, 0],
    bounds=[[0, 0, -6], [1, 10, 6]],
)
print(so_popt)
xv = np.linspace(-3, 6, num=1000)
yv = relu(xv, *so_popt)
ax.plot(xv, yv, color="red", lw=2, linestyle="--", zorder=-10)

sns.scatterplot(
    data=pdf,
    x="avg_enrichment_d2",
    y="sigma_off",
    color="white",
    edgecolor="tab:blue",
    s=40,
    linewidth=2,
    ax=ax,
)

r, p = st.pearsonr(pdf["avg_enrichment_d2"], pdf["sigma_off"])

ax.set_xlim(-3.05, 6.05)
ax.set_xticks([-3, 0, 3, 6])

ax.set_ylim(0, 1)
ax.set_yticks([0, 0.5, 1])

ax.text(x=-2.9, y=0.9, s="Pearson R$^2$={:.2f}".format(r), fontsize=16)

ax.set_xlabel("Act. $\log_2$(ON:OFF)")
ax.set_ylabel(r"$\sigma_{\mathrm{off}}$", rotation="horizontal", ha="right")

fig.savefig("./plots/screen_d2_vs_sigmaoff.pdf", bbox_inches="tight")


def estimate_bprime(d2_score: Union[float, np.ndarray, pd.Series]) -> float:
    return relu(d2_score, *bp_popt)  # type: ignore


def estimate_lambda(d2_score: Union[float, np.ndarray, pd.Series]) -> float:
    return relu(d2_score, *lm_popt)  # type: ignore


def estimate_mu_off(d2_score: Union[float, np.ndarray, pd.Series]) -> float:
    return relu(d2_score, *mu_popt)  # type: ignore


sigma_on_est = np.mean(pdf["sigma_on"])


def estimate_sigma_on(d2_score: Union[float, np.ndarray, pd.Series]) -> float:
    return sigma_on_est


def estimate_sigma_off(d2_score: Union[float, np.ndarray, pd.Series]) -> float:
    return relu(d2_score, *so_popt)  # type: ignore


def estimate_cit_pdf(d2_score: Union[float, np.ndarray, pd.Series]) -> Callable:
    return lambda x: zip_pdf(  # type: ignore
        x,
        estimate_bprime(d2_score),
        estimate_lambda(d2_score),
        estimate_sigma_on(d2_score),
        estimate_sigma_off(d2_score),
        estimate_mu_off(d2_score),
    )


# Let's plot some example histograms

scorevals = [1.5, 2.5, 4.0]


def compute_cit_df(d2_score):
    likelihood = estimate_cit_pdf(d2_score)
    xv = np.linspace(6, 9, num=100)
    yv = likelihood(xv)
    return pd.DataFrame.from_dict(
        {"cit": list(xv), "density": list(yv), "d2_score": [d2_score for x in xv]}
    )


cit_df = pd.concat([compute_cit_df(s) for s in tqdm(scorevals)])

fig, ax = plt.subplots(figsize=(5, 3))

g = sns.lineplot(
    data=cit_df,
    x="cit",
    y="density",
    hue="d2_score",
    palette="flare_r",
    legend="full",
    ax=ax,
)

leg = ax.legend(title="Act. Screen Score", bbox_to_anchor=(1.05, 1.05))

sns.despine(fig)

screen = pd.read_csv("../01_activators_synergy/pairs_baselinesums.csv")

# parameters: bprime, lambda, sigma_off, mu_off
screen["d1_bprime"] = estimate_bprime(screen["d1_med_d2"])
screen["d1_lambda"] = estimate_lambda(screen["d1_med_d2"])
screen["d1_sigma_off"] = estimate_sigma_off(screen["d1_med_d2"])
screen["d1_mu_off"] = estimate_mu_off(screen["d1_med_d2"])

screen["d2_bprime"] = estimate_bprime(screen["d2_med_d2"])
screen["d2_lambda"] = estimate_lambda(screen["d2_med_d2"])
screen["d2_sigma_off"] = estimate_sigma_off(screen["d2_med_d2"])
screen["d2_mu_off"] = estimate_mu_off(screen["d2_med_d2"])

screen["combo_bprime"] = estimate_bprime(screen["avg_enrichment_d2"])
screen["combo_lambda"] = estimate_lambda(screen["avg_enrichment_d2"])
screen["combo_sigma_off"] = estimate_sigma_off(screen["avg_enrichment_d2"])
screen["combo_mu_off"] = estimate_mu_off(screen["avg_enrichment_d2"])

screen["sum_bprime"] = screen["d1_bprime"] + screen["d2_bprime"]
screen["sum_lambda"] = screen["d1_lambda"] + screen["d2_lambda"]
screen["sum_sigma_off"] = screen["d1_sigma_off"] + screen["d2_sigma_off"]
screen["sum_mu_off"] = screen["d1_mu_off"] + screen["d2_mu_off"]

adf = screen[screen["composition"].isin(["A-A", "C-C"])]

fig, ax = plt.subplots(figsize=(3, 2))
sns.scatterplot(
    data=adf,
    x="sum_bprime",
    y="combo_bprime",
    hue="composition",
    linewidth=0,
    marker=".",
    alpha=0.3,
    palette=["#dbc60d", "#999999"],
    ax=ax,
)


combo_less_df = adf[adf["sum_bprime"] >= adf["combo_bprime"]]
combo_more_df = adf[adf["sum_bprime"] < adf["combo_bprime"]]
frac_combo_more = (
    combo_more_df.shape[0] / (combo_more_df.shape[0] + combo_less_df.shape[0]) * 100
)
frac_combo_less = (
    combo_less_df.shape[0] / (combo_more_df.shape[0] + combo_less_df.shape[0]) * 100
)

print("For bprime, frac_combo_more is", frac_combo_more)

ax.text(x=0.45, y=0.85, s="{:.0f}%".format(frac_combo_more))
ax.text(x=1.0, y=0.85, s="{:.0f}%".format(frac_combo_less))


ax.set_xlabel(r"Sum of $\beta'$")
ax.set_xlim(0, 1.5)
ax.set_xticks([0, 0.5, 1.0, 1.5])

ax.set_ylabel(r"Combo $\beta'$")
ax.set_ylim(0, 1.0)
ax.set_yticks([0, 0.5, 1.0])

ax.axline([0, 0], slope=1, color="#777777", linestyle="--", zorder=-10)
leg = ax.legend(
    bbox_to_anchor=(1.0, 1.0),
    fontsize=14,
    markerscale=1.0,
    borderpad=0.3,
    handletextpad=0,
)

plt.savefig("./plots/sumbprime_vs_comboprime.pdf", bbox_inches="tight")

g = sns.JointGrid(
    data=adf,
    x="sum_lambda",
    y="combo_lambda",
    hue="composition",
    palette=["#dbc60d", "#999999"],
    height=4.5,
    marginal_ticks=False,
)

g.plot_joint(
    sns.scatterplot,
    linewidth=0,
    #     marker=".",
    s=40,
    alpha=0.5,
)
# g.plot_marginals(sns.histplot)

sum_less_df = adf[adf["combo_lambda"] > adf["sum_lambda"]]
sum_more_df = adf[adf["combo_lambda"] <= adf["sum_lambda"]]

frac_sum_more = (
    sum_more_df.shape[0] / (sum_more_df.shape[0] + sum_less_df.shape[0]) * 100
)
frac_sum_less = (
    sum_less_df.shape[0] / (sum_more_df.shape[0] + sum_less_df.shape[0]) * 100
)

print("For lambda, frac_sum_more is", frac_sum_more)

g.ax_joint.text(x=1100, y=1425, s="{:.0f}%".format(frac_sum_less))
g.ax_joint.text(x=1325, y=1200, s="{:.0f}%".format(frac_sum_more))

g.ax_joint.axline([0, 0], slope=1, color="#777777", linestyle="--", zorder=-10)

g.ax_joint.set_xlim(-50, 1600)
g.ax_joint.set_ylim(-50, 1600)
g.ax_joint.set_xticks([0, 500, 1000, 1500])
g.ax_joint.set_yticks([0, 500, 1000, 1500])
g.ax_joint.set_xlabel(r"Sum of $\lambda$")
g.ax_joint.set_ylabel(r"Combo $\lambda$")

leg = g.ax_joint.legend(
    markerscale=1.0,
    borderpad=0.2,
    handletextpad=0,
    handlelength=1.5,
    bbox_to_anchor=(0.5, 0.7),
)


plt.savefig("./plots/sumlambda_vs_combolambda.pdf", bbox_inches="tight")

param_titles = [
    r"$\beta'$",
    r"$\lambda$",
    r"$\sigma_\mathrm{Off}$",
    r"$\mu_\mathrm{Off}$",
]

params = ["bprime", "lambda", "sigma_off", "mu_off"]
est_funs = [estimate_bprime, estimate_lambda, estimate_sigma_off, estimate_mu_off]

lower_lims = [0, -50, 0, 6]
upper_lims = [1, 1500, 1, 8]
ticksets = [
    [0, 0.25, 0.5, 0.75, 1.0],
    [0, 500, 1000, 1500],
    [0, 0.25, 0.5, 0.75, 1.0],
    [6, 7, 8],
]

fig, axes = plt.subplots(2, 2, figsize=(8, 8))

for i in tqdm(range(len(params))):
    ax = axes.flat[i]

    xvals = pdf[params[i]]
    yvals = est_funs[i](pdf["avg_enrichment_d2"])
    ax.scatter(
        xvals,
        yvals,
        color="white",
        edgecolor="tab:blue",
        s=40,
        linewidth=2,
    )

    ax.set_title(param_titles[i])

    ax.set_xlim(lower_lims[i], upper_lims[i])
    ax.set_ylim(lower_lims[i], upper_lims[i])

    ax.set_xticks(ticksets[i])
    ax.set_yticks(ticksets[i])

    ax.set_xlabel("Curve Fit")
    ax.set_ylabel("Estimate from Screen")

    r, p = st.pearsonr(xvals, yvals)
    limdelta = upper_lims[i] - lower_lims[i]
    ax.text(
        x=0.05 * limdelta + lower_lims[i],
        y=0.9 * limdelta + lower_lims[i],
        s="Pearson R$^2$={:.2f}".format(r),
        fontsize=16,
    )

plt.tight_layout()

plt.savefig("./plots/screen_params_vs_fit_values.pdf", bbox_inches="tight")


# from params, compute fraction on, compare to avg_enrichment_d2???
def compute_fon_from_params(bprime, lmbda, sigma_on, sigma_off, mu_off):
    def likelihood(c):
        return lambda c: zip_pdf(c, bprime, lmbda, sigma_on, sigma_off, mu_off)

    f_on = 1 - (
        scipy.integrate.quad(likelihood, 6, 7)[0]  # type: ignore
        / scipy.integrate.quad(likelihood, 6, 9)[0]  # type: ignore
    )
    return f_on


def compute_fon_from_screen(d2_score):
    likelihood = estimate_cit_pdf(d2_score)
    f_on = 1 - (
        scipy.integrate.quad(likelihood, 6, 7)[0]  # type: ignore
        / scipy.integrate.quad(likelihood, 6, 9)[0]  # type: ignore
    )
    return f_on


# Beta vs fraction on

scores = np.linspace(-3, 6, num=100)
betas = [estimate_bprime(s) for s in tqdm(scores)]
lmbdas = [estimate_lambda(s) for s in tqdm(scores)]
fons = [compute_fon_from_screen(s) for s in tqdm(scores)]

sbfdf = pd.DataFrame.from_dict(
    {"score": scores, "beta": betas, "lambda": lmbdas, "fon": fons}
)


fig, ax = plt.subplots(1, 2, figsize=(10, 4))
sns.scatterplot(
    data=sbfdf,
    x="beta",
    y="fon",
    hue="score",
    palette="flare",
    legend="brief",
    ax=ax[0],  # type: ignore
)

sns.scatterplot(
    data=sbfdf,
    x="lambda",
    y="fon",
    hue="score",
    palette="flare",
    legend="brief",
    ax=ax[1],  # type: ignore
)

ax[0].legend().set_visible(False)  # type: ignore
leg = ax[1].legend(title="Act. Screen Score", bbox_to_anchor=(1.1, 1.05), markerscale=1.0)  # type: ignore

ax[0].set_xlabel(r"$\beta'$")  # type: ignore
ax[0].set_ylabel("Est. Fraction On")  # type: ignore
ax[0].set_xlim(0, 1)  # type: ignore
ax[0].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])  # type: ignore

ax[1].set_xlabel(r"$\lambda$")  # type: ignore
ax[1].set_ylabel("")  # type: ignore
ax[1].set_xlim(-50, 1500)  # type: ignore
ax[1].set_xticks([0, 500, 1000, 1500])  # type: ignore

for a in ax:  # type: ignore
    a.set_ylim(0, 1)
    a.set_yticks([0, 0.5, 1])

plt.tight_layout()


ax[0].axline((0, 0), slope=1, color="#777777", linestyle="--", lw=2, zorder=-10)  # type: ignore
ax[1].axline((0, 0), slope=1 / 1500, color="#777777", linestyle="--", lw=2, zorder=-10)  # type: ignore

plt.savefig("./plots/est_fon_from_bprime_lamba.pdf", bbox_inches="tight")
