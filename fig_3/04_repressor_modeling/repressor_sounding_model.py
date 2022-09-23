from typing import List
from autograd.scipy.stats import norm
from autograd import value_and_grad
from scipy.optimize import minimize
from scipy.stats import norm as norm_
import numpy as np


def compute_fractions(xdata: np.ndarray) -> List:  # type: ignore
    times = xdata[:, 0]
    log_cits = xdata[:, 1]


def negative_log_likelihood(params, obs):  # assume 3 peaks
    mu1, ls1, w1, mu2, ls2, w2, mu3, ls3, w3 = params
    s1, s2, s3 = np.exp([ls1, ls2, ls3])
    return -np.log(
        w1 * norm.pdf(obs, mu1, s1)
        + w2 * norm.pdf(obs, mu2, s2)
        + w3 * norm.pdf(obs, mu3, s3)
    )

# generate some fake data. Success will be if we recover these parameters. 
obs = np.r_[
    norm_(loc=6.5, scale=0.25).rvs(200),
    norm_(loc=7.5, scale=0.5).rvs(500),
    norm_(loc=8.5, scale=1).rvs(300)
]


results = minimize(
    value_and_grad(negative_log_likelihood), # see autograd docs.
    x0 = np.array([6, 0, 0.33, 7, 0, 0.33, 8, 0, 0.33]), # initial value
    args=(obs,),
    jac=True,
    bounds=(
        (6, 7),    # mu1 (you mentioned the constraints on the means)
        (None, None), # log_sigma1 is unbounded
        (0, 1), # w1
        (7, 8),    # mu2 (you mentioned the constraints on the means)
        (None, None), # log_sigma2 is unbounded
        (0, 1),
        (8, 9),
        (None, None),
        (0, 1)        # the weight param should be between 0 and 1
    )
)

print(results.x)

