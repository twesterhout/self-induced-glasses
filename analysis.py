import glob
import math
import os
import re

import numpy as np
from numpy.typing import NDArray, ArrayLike
import scipy


def autocorr_function(x: NDArray[np.float64]) -> NDArray[np.float64]:
    r"""Estimate the normalized autocorrelation function of a 1D array."""
    if x.ndim != 1:
        raise ValueError("x has wrong shape: {}; expected a 1D array".format(x.shape))
    n = 1 << math.ceil(math.log2(len(x)))
    f = np.fft.fft(x - np.mean(x), n=2 * n)
    autocorr = np.fft.ifft(f * np.conj(f))[: len(x)].real
    autocorr /= autocorr[0]
    return autocorr


def extract_temperature(filename: str) -> float:
    return float(re.match(r".*measurements_T=(.+)\.csv", filename).group(1))


def get_temperature_set(prefix):
    return set(
        extract_temperature(filename)
        for filename in glob.glob("{}/measurements_T=*.csv".format(prefix))
    )

def load_measurements_with_disorder(prefix="ising-4/*"):
    temperatures = get_temperature_set(prefix)
    for T in sorted(list(temperatures)):
        data = np.vstack(
            [
                np.loadtxt(filename, delimiter=",")
                for filename in glob.glob("{}/measurements_T={:.4f}.csv".format(prefix, T))
            ]
        )
        assert data.shape[1] == 7
        yield (T, data)


def compute_binder_ratio(prefix: str = "ising-4/*", output: str = "ising-4"):
    os.makedirs(output, exist_ok=True)

    def binder_fn(data):
        q = data[:, 6]
        return 0.5 * (3 - np.mean(q ** 4) / np.mean(q ** 2) ** 2)

    results = np.asarray(
        [(T, binder_fn(data)) for T, data in load_measurements_with_disorder(prefix)]
    )
    np.savetxt(
        "{}/binder_ratio.csv".format(output),
        results,
        delimiter=",",
    )

def compute_fake_susceptibility(n: int, prefix: str = "ising-4/*", output: str = "ising-4"):
    os.makedirs(output, exist_ok=True)

    def susceptibility_fn(data):
        q = data[:, 6]
        return n * np.var(np.abs(q))

    results = np.asarray(
        [(T, susceptibility_fn(data)) for T, data in load_measurements_with_disorder(prefix)]
    )
    np.savetxt(
        "{}/fake_susceptibility.csv".format(output),
        results,
        delimiter=",",
    )

def compute_susceptibility_scaling(n: int, prefix: str = "ising-4/*", output: str = "ising-4"):
    os.makedirs(output, exist_ok=True)

    def susceptibility_fn(data):
        q = data[:, 6]
        return n * np.var(np.abs(q))

    results = np.asarray(
        [(T, susceptibility_fn(data)) for T, data in load_measurements_with_disorder(prefix)]
    )
    np.savetxt(
        "{}/fake_susceptibility.csv".format(output),
        results,
        delimiter=",",
    )


def compute_density_with_disorder(
    prefix: str = "ising-4/*", output: str = "ising-4/P_q", eps: float = 5e-2, points: int = 500
):
    os.makedirs(output, exist_ok=True)
    xs = np.linspace(-1, 1, points)
    for T, data in load_measurements_with_disorder(prefix):
        kernel = scipy.stats.gaussian_kde(data[:, 6], bw_method=eps)
        np.savetxt(
            "{}/P_q_T={:.4f}.csv".format(output, T),
            np.vstack([xs, kernel(xs)]).T,
            delimiter=",",
        )


def compute_densities(eps: float):
    xs = np.linspace(-1, 1, 500)
    for filename in glob.glob("measurements_T=*.csv"):
        T = extract_temperature(filename)
        table = np.loadtxt(filename, delimiter=",")
        kernel = scipy.stats.gaussian_kde(table[:, 6], bw_method=eps)
        np.savetxt("P_q_T={:.4f}.csv".format(T), np.vstack([xs, kernel(xs)]).T, delimiter=",")


def main():
    for filename in glob.glob("measurements_T=*.csv"):
        T = extract_temperature(filename)
        table = np.loadtxt(filename, delimiter=",")
        f1 = autocorr_function(table[:, 0])
        f2 = autocorr_function(table[:, 1])
        np.savetxt("autocorr_T={:.4f}.csv".format(T), np.vstack([f1, f2]).T, delimiter=",")


if __name__ == "__main__":
    main()
