import os
from typing import List, Optional, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.integrate
import scipy.optimize
import scipy.stats as st
import seaborn as sns
from matplotlib.lines import Line2D
from sklearn.mixture import GaussianMixture
from tqdm import tqdm

mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["font.family"] = "sans-serif"
plt.rcParams["figure.figsize"] = (5, 3)  # type: ignore
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
sns.set_context("talk", font_scale=1)

DEBUG = False
min_cit = 5.5
max_cit = 9.5
max_mrna = 600 * (max_cit - min_cit)


def printdb(s: str) -> None:
    if DEBUG:
        print(s)


def calculate_a(t: np.ndarray, ks: float, tlag: float) -> np.ndarray:
    """
    calculate_a computes the active fraction using the 3-state model

    Args:
        t (np.ndarray): time since recruitment began
        ks (float): rate of reversible silencing
        tlag (float): lag time before silencing

    Returns:
        np.ndarray: array of fractions on between 0 and 1 of active cells
    """
    return 1.0 * (t < tlag) + np.exp(-1 * ks * (t - tlag)) * (t >= tlag)


def prob_cit_on_simple(
    log_cits: np.ndarray,
    times: np.ndarray,
    tlag: float,
    lmbda: float,
    sigma_on: float,
    gamma: float,
    mu: float,
    sigma_off: float,
) -> np.ndarray:
    """
    prob_cit_on_simple computes the pdf of citrine in the ON population
    (either active or silent)

    Args:
        log_cits (np.ndarray): log 10 of citrine levels
        times (np.ndarray): time in days corresponding to citrine levels
        tlag (float): lag time before silencing begins
        lmbda (float): net ratio of mRNA production to degradation
        sigma_on (float): standard deviation of the ON population
        gamma (float): rate of mRNA degradation + dilution
        mu (float): mean of OFF population
        sigma_off (float): standard deviation of OFF population

    Returns:
        np.ndarray: array of probability densities corresponding to citrines
    """
    mrna_levels = np.arange(0, 600 * (max_cit - mu))
    mrna_means = lmbda * (times < tlag) + lmbda * np.exp(-gamma * (times - tlag)) * (
        times >= tlag
    )
    sigmas = np.array([sigma_on if t <= 2 else sigma_off for t in times])
    p_m = st.poisson(mrna_means).pmf(mrna_levels.reshape(-1, 1))
    p_c = st.norm(((mrna_levels.reshape(-1, 1) / 600) + mu).reshape(-1, 1), sigmas).pdf(
        log_cits
    )
    res = np.sum(p_m * p_c, axis=0)
    return res  # asumes unif distribution over mrna levels


def telegraph_3sm_pdf(
    xdata: np.ndarray,
    ks: float,
    tlag: float,
    bprime: float,
    lmbda: float,
    s_on: float,
    gamma: float,
    u_off: float,
    s_off: float,
) -> np.ndarray:
    """
    telegraph_3sm_pdf computes the pdf using the 3-state model

    Args:
        xdata (np.ndarray): Mx2 array of time being the 1st col and cit being the 2nd
        ks (float): rate of silencing
        tlag (float): lag time before silencing begins
        bprime (float): fraction of cells ON at equilibrium
        lmbda (float): ratio of mRNA production to degradation
        s_on (float): standard deviation of ON population
        gamma (float): net rate of mRNA degradation + dilution
        u_off (float): mean of OFF population
        s_off (float): standard deviation of OFF population

    Returns:
        np.ndarray: probability density of citrine, time tuples
    """
    times = xdata[:, 0]
    log_cits = xdata[:, 1]

    # cit_range = (np.arange(0, max_mrna) / 600) + min_cit
    # est_mrnas = np.rint(u_off + ((log_cits - u_off) * 600))

    p_active = bprime * calculate_a(times, ks, tlag)
    p_silent = bprime * (1 - p_active)
    p_off = 1 - bprime

    p_cit_active = prob_cit_on_simple(
        log_cits, np.tile(0, len(log_cits)), tlag, lmbda, s_on, gamma, u_off, s_off
    )
    if ks != 0:
        # p_cit_silent = prob_cit_on_silent(
        #     log_cits, times, ks, lmbda, s_on, gamma, u_off, s_off
        # )
        p_cit_silent = prob_cit_on_simple(
            log_cits, times, tlag, lmbda, s_on, gamma, u_off, s_off
        )
    else:
        p_cit_silent = np.array([0 for _ in log_cits])
    p_cit_off = st.norm(u_off, s_off).pdf(log_cits)

    return p_cit_active * p_active + p_cit_silent * p_silent + p_off * p_cit_off


def get_ks_tlag_gamma(
    times: np.ndarray, log_cits: np.ndarray, u_off: float
) -> List[float]:
    """
    get_ks_tlag_gamma computes ks, tlag, and gamma from citrine data

    Args:
        times (np.ndarray): list of times in days
        log_cits (np.ndarray): list of log citrine values
        u_off (float): mean of OFF population

    Returns:
        List[float]: list of [ks, tlag, gamma]
    """
    dq = pd.DataFrame.from_dict({"time": times, "citrine": log_cits})

    # get fractions on
    dq["on"] = dq["citrine"] >= 8
    dx = dq.sort_values(by="time").groupby("time").mean().reset_index()
    ont = list(dx["time"])
    onv = list(dx["on"])

    def get_frac_on(cits: np.ndarray, d: float) -> float:
        """
        get_frac_on computes the fraction of cells on

        Args:
            cits (np.ndarray): list of citrine values
            d (float): day

        Returns:
            float: fraction of cells on
        """
        # i do not care that this code is repeated at the moment
        gm = GaussianMixture(n_components=2).fit(cits.reshape(-1, 1))
        m1, m2 = gm.means_
        w1, w2 = gm.weights_
        if m1 < 8 and m2 < 8:  # everything is "off"
            return onv[ont.index(d)]
        else:  # return frac in on peak
            return w1 if m1 > m2 else w2

    def get_loc_on(cits: np.ndarray, d) -> float:
        """
        get_loc_on returns the location of the ON peak

        Args:
            cits (np.ndarray): list of citrine values
            d (_type_): day

        Returns:
            float: ON peak location
        """
        gm = GaussianMixture(n_components=2).fit(cits.reshape(-1, 1))
        m1, m2 = gm.means_
        w1, w2 = gm.weights_
        m1 = m1[0]
        m2 = m2[0]
        m_min = min(m1, m2)
        m_max = max(m1, m2)
        w_min = w1 if m_min == m1 else w2
        w_max = w1 if m_max == m1 else w2
        if w_max < 0.1:
            return m_min
        elif w_min < 0.1:
            return m_max
        else:
            return m_max

    tvals = sorted(list(set(list(dq["time"]))))
    fvals = [get_frac_on(np.array(dq[dq["time"] == d]["citrine"]), d) for d in tvals]
    cvals = [get_loc_on(np.array(dq[dq["time"] == d]["citrine"]), d) for d in tvals]

    start_fon = fvals[0]

    def fraction_off_curve(t, ks, tlag):
        return start_fon * (1 * (t < tlag) + np.exp(-ks * (t - tlag)) * (t >= tlag))

    def peak_decay_curve(t, gamma, tlag):
        return np.max(cvals) * (t < tlag) + (
            u_off + (np.max(cvals) - u_off) * np.exp(-gamma * (t - tlag))
        ) * (t >= tlag)

    fopt, _ = scipy.optimize.curve_fit(
        fraction_off_curve,
        tvals,
        fvals,
        p0=[5, 1],
        bounds=[[0, 0], [15, 2]],
        ftol=1e-14,
    )
    ks, t1 = fopt

    copt, _ = scipy.optimize.curve_fit(
        lambda x, p: peak_decay_curve(x, p, t1),
        tvals,
        cvals,
        p0=[0.5],
        bounds=[[0.4], [1.0]],
    )
    g = copt[0]

    return [ks, t1, g]


def get_fit_params(
    times: np.ndarray, log_cits: np.ndarray, plasmid: float
) -> List[float]:
    """
    get_fit_params returns optimal fit parameters for the telegraph model given data

    Args:
        times (np.ndarray): numpy array of times
        log_cits (np.ndarray): numpy array of log citrine measurements
        plasmid (float): number corresponding to the plasmid of interest

    Returns:
        List[float]: list of optimal parameters, in the order:
            ks, tlag, bprime, lambda, sigma_on, gamma, mu_off, sigma_off
    """
    xdata = np.transpose(np.array([times, log_cits]))

    # first, fit beta and lambda to day 0 data
    zxdata = xdata[xdata[:, 0] == 0]
    zkde = st.gaussian_kde(dataset=zxdata[:, 1])
    zydata = zkde.evaluate(zxdata[:, 1])

    def zero_fn(x, bp, lmbda, s_on):
        return telegraph_3sm_pdf(x, 0, 10, bp, lmbda, s_on, 0, 6.4, 0.25)

    # printdb("Fitting beta, lambda")
    zopt, _ = scipy.optimize.curve_fit(
        f=zero_fn,
        xdata=zxdata,
        ydata=zydata,
        p0=[0.5, 1200, 0.5],
        bounds=[[0, 600, 0], [1, max_mrna, 2]],
    )
    fit_b, fit_l, fit_son = zopt
    printdb(
        "\tFound beta={:.2f}, lambda={:.2f}, sigma_on={:.2f}".format(
            fit_b, fit_l, fit_son
        )
    )

    # next, fit mu and sigma to last day data less than 1e7
    lxdata = xdata[xdata[:, 0] == np.max(xdata[:, 0])]
    lxdata = lxdata[lxdata[:, 1] <= 7]
    lkde = st.gaussian_kde(dataset=lxdata[:, 1])
    lydata = lkde.evaluate(lxdata[:, 1])

    def off_fn(x, mu, sigma):
        # here we pick super strong silencing params
        return telegraph_3sm_pdf(x, 10, 0, 0, max_mrna, 1e-200, 10, mu, sigma)

    # printdb("Fitting mu, sigma")
    lopt, _ = scipy.optimize.curve_fit(
        f=off_fn,
        xdata=lxdata,
        ydata=lydata,
        p0=[6.3, 0.5],
        bounds=[[5.5, 0.1], [7, 1.0]],
    )
    fit_m, fit_soff = lopt
    printdb("\tFound mu={:.2f} and sigma_off={:.2f}".format(fit_m, fit_soff))

    # last, fit ks, tlag, and gamma
    # printdb("Fitting ks, tlag, gamma")
    fit_k, fit_t, fit_g = get_ks_tlag_gamma(times, log_cits, fit_m)
    fit_k = fit_k if fit_k > 0.1 else 0

    printdb(
        "\tFound ks={:.2f}, tlag={:.2f}, and gamma={:.2f}".format(fit_k, fit_t, fit_g)
    )

    return [fit_k, fit_t, fit_b, fit_l, fit_son, fit_g, fit_m, fit_soff]


def import_repression_data() -> pd.DataFrame:
    """
    import_repression_data imports repression data from the two measurements

    Returns:
        pd.DataFrame: DF of all cells measured, gated for mch+ and singlets
    """
    rep_round1 = pd.read_csv("./data/rep_rd1_all_cells_live_mch_gated.csv")
    rep_round1 = rep_round1[rep_round1["treatment"] == "none"]
    rep_round1["date"] = "2021.11.17"

    rep_round2 = pd.read_csv("./data/repadd_stapl_all_cells_mch_live.csv")
    rep_round2 = rep_round2[rep_round2["asv"] == 0]
    rep_round2 = rep_round2.drop("asv", axis=1)
    rep_round2["date"] = "2022.04.22"

    rep_df = pd.concat([rep_round1, rep_round2])
    rep_df = rep_df.replace([np.inf, -np.inf], np.nan)
    rep_df = rep_df.dropna(subset=["mCitrine-A"])
    return rep_df


def get_test_histogram_data(
    plasmid: int = 217, df: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """
    get_test_histogram_data gets test histogram data for a given plasmid

    Args:
        plasmid (int, optional): plasmid number, defaults to 217.

    Returns:
        pd.DataFrame: dataframe of all +dox cells for that plasmid
    """
    if df is None:
        df = import_repression_data()
    df = df[df["plasmid"] == plasmid]
    df = df[df["dox"] == 1000]
    return df


def plot_param_breakdown(params: List[float]) -> plt.Axes:
    """
    plot_param_breakdown Plots a breakdown of ON/ACTIVE/SILENT states over time

    Args:
        params (List[float]): list of parameters:
            ks, tlag, bprime, lambda, son, gamma, u, soff

    Returns:
        plt.Axes: axis object of plot
    """
    ks, tlag, bprime, lmbda, son, gamma, u, soff = params

    cits = np.linspace(5.5, 9.5, num=1000)

    def get_times(day):
        return np.array([day for _ in cits])

    def p_active(day):
        return calculate_a(get_times(day), ks, tlag)

    def p_silent(day):
        return 1 - p_active(day)

    p_off = 1 - bprime

    def get_off(day):
        return p_off * st.norm(u, soff).pdf(cits)

    def get_silent(day):
        return (
            p_silent(day)
            * bprime
            * prob_cit_on_simple(cits, get_times(day), tlag, lmbda, son, gamma, u, soff)
            # * prob_cit_on_silent(cits, get_times(day), ks, lmbda, son, gamma, u, soff)
        )

    def get_active(day):
        return (
            p_active(day)
            * bprime
            * prob_cit_on_simple(cits, get_times(0), tlag, lmbda, son, gamma, u, soff)
        )

    def get_total(day):
        return get_active(day) + get_silent(day) + get_off(day)

    daylist = [0, 1, 2, 3, 4, 5]
    nrows = int(max(1, len(daylist) / 3))
    ncols = 3
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5 * ncols, 3 * nrows))
    for i, a in tqdm(enumerate(ax.flat), total=len(ax.flat)):
        day = daylist[i]

        off = get_off(day)
        sil = get_silent(day)
        act = get_active(day)
        tot = get_total(day)

        xvals = np.tile(cits, 4)
        yvals = np.concatenate([off, sil, act, tot])
        labs = np.concatenate(
            [
                np.array(["off" for _ in off]),
                np.array(["sil" for _ in sil]),
                np.array(["act" for _ in act]),
                np.array(["tot" for _ in tot]),
            ]
        )

        pdf = pd.DataFrame.from_dict({"cit": xvals, "prob": yvals, "type": labs})
        # pdf = pdf[pdf["type"].isin(["off", "sil"])]

        sns.lineplot(data=pdf, x="cit", y="prob", hue="type", palette="Set2", ax=a)

        a.set_xlim(5.5, 9.5)
        a.set_xticks([6, 7, 8, 9])
        a.set_yticks([])

    return ax


def fit_and_plot(
    plasmid_df: Union[int, Optional[pd.DataFrame]] = None,
    params: Optional[List[float]] = None,
):
    """
    fit_and_plot fits the telegraph model to plasmid data and plots the fit

    Args:
        plasmid_df (Union[int, Optional[pd.DataFrame]], optional): plasmid dataframe.
                Defaults to None.
        params (Optional[List[float]], optional): list of parameters. Defaults to None.

    Returns:
        [plt.Figure, pd.DataFrame]: dataframe of parameter fits, plot of the fit
    """
    if plasmid_df is None:
        printdb("Reading in data with default plasmid")
        plasmid_df = get_test_histogram_data()
    elif isinstance(plasmid_df, int):
        printdb("Reading in data for plasmid " + str(plasmid_df))  # type: ignore
        plasmid_df = get_test_histogram_data(
            plasmid_df
        )  # enabling some real abuse here
    else:
        printdb("User provided dataframe, no need to load")
    plasmid = list(plasmid_df["plasmid"])[0]
    descr = list(plasmid_df["description"])[0]

    if params is None:
        printdb("Getting fit parameters")
        plasmid_df = plasmid_df[plasmid_df["mCitrine-A"] > 0]
        times = np.array(plasmid_df["day"])
        lcits = np.log10(np.array(plasmid_df["mCitrine-A"]))
        popt = get_fit_params(times, lcits, plasmid)
    else:
        printdb("Given parameters, going directly to plotting")
        popt = params

    daylist = sorted(list(set(list(plasmid_df["day"]))))
    nrows = int(max(1, len(daylist) / 3))
    ncols = 3
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5 * ncols, 3 * nrows))
    for i, a in enumerate(ax.flat):
        day = int(daylist[i])
        # print(f"Working on day {day:d}")

        def popt_fun(log_cits: np.ndarray) -> np.ndarray:
            days = np.array([day for _ in log_cits])
            xdata = np.transpose(np.array([days, log_cits]))
            return telegraph_3sm_pdf(xdata, *popt)

        sns.kdeplot(
            data=plasmid_df[plasmid_df["day"] == day],
            x="mCitrine-A",
            color="tab:red",
            log_scale=True,
            ax=a,
        )

        log_cits = np.linspace(min_cit, max_cit)
        yvals = popt_fun(log_cits)
        xvals = np.power(10, log_cits)
        a.plot(xvals, yvals, color="tab:blue", linestyle="--")

        a.set_title("Day " + str(int(day)))
        a.set_xscale("log")
        a.set_xlim(3.16e5, 3.16e9)
        a.set_xticks([1e6, 1e7, 1e8, 1e9])

    custom_lines = [
        Line2D([0], [0], color="tab:red", lw=4),
        Line2D([0], [0], color="tab:blue", lw=4),
    ]
    fig.axes[0].legend(
        custom_lines,
        ["Data", "Fit"],
        loc="upper left",
        #        bbox_to_anchor=(1.15, -7),
    )

    sns.despine(fig)

    plt.tight_layout()
    # plt.show()

    fig.suptitle(descr)

    ks, tl, bp, lm, son, y, uo, soff = popt
    param_df = pd.DataFrame.from_dict(
        {
            "plasmid": [list(plasmid_df["plasmid"])[0]],
            "ks": [ks],
            "tlag": [tl],
            "bprime": [bp],
            "lambda": [lm],
            "s_on": [son],
            "gamma": [y],
            "u_off": [uo],
            "s_off": [soff],
        }
    )
    return [fig, param_df]


def fit_all() -> pd.DataFrame:
    df = import_repression_data()
    plasmid_list = sorted(list(set(list(df["plasmid"]))))

    def get_descr(p):
        return list(df[df["plasmid"] == p]["description"])[0]

    descr_list = list(map(get_descr, plasmid_list))

    print("Running telegraph model fits for", len(plasmid_list), "plasmids")

    dfs = []
    for i, p in tqdm(enumerate(plasmid_list), total=len(plasmid_list)):
        descr = descr_list[i]
        full_p_data = get_test_histogram_data(p, df)
        sample_size = np.min(list(full_p_data.groupby("day").count()["dox"])) - 1
        p_data = full_p_data.groupby("day").sample(n=sample_size)
        fig, pdf = fit_and_plot(p_data)
        fig.savefig("./plot_telegraph/" + descr + "_fit.pdf", bbox_inches="tight")
        dfs.append(pdf)
        plt.close(fig)

    param_df = pd.DataFrame(pd.concat(dfs))
    param_df.to_csv("./data/telegraph_parameters.csv")

    return param_df


def make_ks_comp_plot(param_df: pd.DataFrame):
    """
    make_ks_comp_plot Makes a plot comparing the telegraph model ks to the prior
    validation fits

    Args:
        param_df (pd.DataFrame): dataframe of ks and prior validation ks parameters

    Returns:
        mpl.Figure: makes plots, returns a Figure object
    """
    fig, ax = plt.subplots(figsize=(3, 3))

    plot_df = param_df.dropna(subset=["ks", "ks_validation"])

    _ = sns.regplot(
        data=plot_df,
        x="ks",
        y="ks_validation",
        scatter=False,
        line_kws={"color": "tab:red", "linestyle": "--", "zorder": -10},
    )

    _ = sns.scatterplot(
        data=plot_df,
        x="ks",
        y="ks_validation",
        color="white",
        edgecolor="tab:blue",
        s=40,
        linewidth=2,
        ax=ax,
    )

    r, p = st.pearsonr(plot_df["ks"], plot_df["ks_validation"])
    ax.text(0.1, 3.5, "Pearson\n$R={:.2f}$".format(r))
    ax.set_xlim(0, 4.5)
    ax.set_xticks([0, 2, 4])
    ax.set_xlabel("Telegraph Model $k_s$")
    ax.set_ylim(0, 4.5)
    ax.set_yticks([0, 2, 4])
    ax.set_ylabel("3-State Model $k_s$")

    return fig


def get_param_df() -> List[pd.DataFrame]:
    """
    get_param_df gets dataframe of validation parameter fits with telegraph model

    Note that if the file at `./data/telegraph_parameters.csv` exists, it will
    simply load that. To force re-calculation, delete the file.

    Returns:
        List[pd.DataFrame]:
            [Dataframe of parameter estimates, Dataframe of prior validation data]
    """
    if not os.path.exists("./data/telegraph_parameters.csv"):
        param_df = fit_all()
    param_df = pd.read_csv("./data/telegraph_parameters.csv")
    val_1_p_df = pd.read_csv("./data/paramed_concats.csv")
    val_2_p_df = pd.read_csv("./data/parameter_df_r372.csv")

    val_df = pd.concat([val_1_p_df, val_2_p_df])
    # val_df["validation_ks"] = val_df["ks"]
    # val_df["validation_tlag"] = val_df["tlag"]
    # val_df = val_df[["plasmid", "validation_ks", "validation_tlag"]]

    param_df = (
        param_df.set_index("plasmid")
        .join(val_df.set_index("plasmid"), how="left", rsuffix="_validation")
        .reset_index()
    )
    return [param_df, val_df]


def make_ks_screen_scatter(joined_df):
    """
    make_ks_screen_scatter Builds plot comparing ks to screen score

    Args:
        joined_df (_type_): dataframe containing ks and screen score data

    Returns:
        mpl.Figure: figure object of plot
    """
    fig, ax = plt.subplots(figsize=(3, 3))

    _ = sns.regplot(
        data=joined_df,
        x="avg_enrichment_d5",
        y="ks",
        scatter=False,
        line_kws={"color": "tab:red", "linestyle": "--", "zorder": -10},
    )

    _ = sns.scatterplot(
        data=joined_df,
        x="avg_enrichment_d5",
        y="ks",
        color="white",
        edgecolor="tab:blue",
        s=40,
        linewidth=2,
        ax=ax,
    )

    r, p = st.pearsonr(joined_df["avg_enrichment_d5"], joined_df["ks"])
    ax.text(-2.5, 3.35, "Pearson\nR$={:.2f}$".format(r))
    ax.set_xlim(-5.5, 2.0)
    ax.set_xticks([-4, -2, 0, 2])
    ax.set_xlabel("Screen $\log_2$(ON:OFF)") 
    ax.set_ylim(-0.2, 4.5)
    ax.set_yticks([0, 2, 4])
    ax.set_ylabel("Telegraph Model $k_s$")

    fig.savefig("./plot_telegraph/ks_screen_comparison.pdf", bbox_inches="tight")

    return fig


if __name__ == "__main__":
    param_df, val_df = get_param_df()
    fig = make_ks_comp_plot(param_df)
    fig.savefig("./plot_telegraph/ks_comparison.pdf", bbox_inches="tight")
    plt.close(fig)

    screen_df = pd.read_csv(
        "../../../github_repo/fig_3/"
        + "01_repressor_additivity/pairs_baselinesums.csv"
    )
    d1g = list(screen_df["d1_Gene"])
    d2g = list(screen_df["d2_Gene"])
    dpairs = [str(d1) + " - " + str(d2) for d1, d2 in zip(d1g, d2g)]
    screen_df["description"] = dpairs

    joined_df = (
        param_df.set_index("description")
        .join(screen_df.set_index("description"), how="left", rsuffix="screen")
        .reset_index()
    )
    joined_df = joined_df.dropna(
        subset=["description", "avg_enrichment_d5", "ks_validation"]
    )
    joined_df = joined_df[~joined_df["composition"].str.contains("C")]
    f = make_ks_screen_scatter(joined_df)
    plt.close(f)
