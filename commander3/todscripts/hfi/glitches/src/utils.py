import numpy as np


def chi2(res):
    sigma0 = np.std(res)
    chi2 = np.sum((res / sigma0) ** 2)
    return chi2

def normalized_chi2(res):
    sigma0 = np.std(res)
    chi2 = np.sum((res / sigma0) ** 2)
    N = len(res)
    norm_chi2 = (chi2 - N) / np.sqrt(2 * N)
    return norm_chi2