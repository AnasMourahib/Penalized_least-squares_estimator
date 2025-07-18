import numpy as np
#import clef_algo as clf
import itertools as it


def extreme_points_bin(x_raw, k=None, R=None, eps=None, without_zeros=None):
    """
        Output:
            -Binary matrix : if k -> kth largest points on each column
                             if R -> points above R
        (Same output if columns are normalized to Pareto(1) distribution
         and R=n_sample/(k + 1))
    """
    n_sample, n_dim = np.shape(x_raw)
    if k:
        mat_rank = np.argsort(x_raw, axis=0)[::-1]
        x_bin = np.zeros((n_sample, n_dim))
        for j in range(n_dim):
            x_bin[mat_rank[:k, j], j] = 1
    if R and not eps:
        x_bin = 1.*(x_raw > R)
    if eps:
        x_bin = 1.*(x_raw[np.max(x_raw, axis=1) > R] > R*eps)
    if without_zeros:
        x_bin = x_bin[np.sum(x_bin, axis=1) > 0]

    return x_bin


def rank_transformation(x_raw):
    n_sample, n_dim = np.shape(x_raw)
    mat_rank = np.argsort(x_raw, axis=0)[::-1]
    x_rank = np.zeros((n_sample, n_dim))
    for i in range(n_dim):
        x_rank[mat_rank[:, i], i] = np.arange(n_sample) + 1
    x_pareto = n_sample/x_rank

    return x_pareto





def r(x_bin, alpha, k):
    """Dependence function


    Parameters
    ----------
    x_bin : ndarray, shape (n_sample, n_dim)
    alpha : list
        List of features
    k : int

    Returns
    -------
    r : float


    """
    return np.sum(np.sum(x_bin[:, alpha], axis=1) == len(alpha))/float(k)


def init_freq(x_bin_k, k, f_min):
    n_dim = np.shape(x_bin_k)[1]
    alphas = []
    for (i, j) in it.combinations(range(n_dim), 2):
        alpha = [i, j]
        r_alph = r(x_bin_k, alpha, k)
        if r_alph > f_min:
            alphas.append(alpha)

    return alphas







def rhos_alpha_pairs(x_bin, alpha, k):
    """Computes rho(i,j) for each (i,j) in alpha.


    Parameters
    ----------
    x_bin : ndarray, shape (n_sample, n_dim)
    alpha : list
        List of features
    k : int

    Returns
    -------
    rhos_alpha : dict


    """
    rhos_alpha = {}
    for (i, j) in it.combinations(alpha, 2):
        rhos_alpha[i, j] = r(x_bin, [i, j], k)

    return rhos_alpha


def partial_matrix(x_bin_base, x_bin_partial, j):
    """Returns x_bin_base with its jth colomn replaced by the jth column of
    x_bin_partial.


    Parameters
    ----------
    x_bin_base : ndarray, shape (n_sample, n_dim)
    x_bin_partial : ndarray, shape (n_sample, n_dim)
    j : int

    Returns
    -------
    x_bin_copy : ndarray, shape (n_sample, n_dim)


    """
    x_bin_copy = np.copy(x_bin_base)
    x_bin_copy[:, j] = x_bin_partial[:, j]

    return x_bin_copy


def r_partial_derv_centered(x_bin_k, x_bin_kp, x_bin_km, alpha, k):
    """
        Output:
            - dictionary : {j: derivative of r in j}
    """
    r_p = {}
    for j in alpha:
        x_r = partial_matrix(x_bin_k, x_bin_kp, j)
        x_l = partial_matrix(x_bin_k, x_bin_km, j)
        r_p[j] = 0.5*k**0.25*(r(x_r, alpha, k) - r(x_l, alpha, k))

    return r_p
