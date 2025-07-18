import numpy as np


#############
# Functions #
#############


def damex_0(x_bin):
    """
    binary x_bin -> nb of points per subfaces
    """
    n_sample, n_dim = np.shape(x_bin)
    n_extr_feats = np.sum(x_bin, axis=1)
    n_shared_feats = np.dot(x_bin, x_bin.T)
    exact_extr_feats = (n_shared_feats == n_extr_feats) * (
        n_shared_feats.T == n_extr_feats).T
    feat_non_covered = set(range(n_sample))
    samples_nb = {}
    for i in range(n_sample):
        feats = list(np.nonzero(exact_extr_feats[i, :])[0])
        if i in feat_non_covered:
            feat_non_covered -= set(feats)
            if n_extr_feats[i] > 1:
                samples_nb[i] = len(feats)
    ind_sort = np.argsort(list(samples_nb.values()))[::-1]
    alphas = [list(np.nonzero(x_bin[list(samples_nb.keys())[i], :])[0])
              for i in ind_sort]
    mass = [list(samples_nb.values())[i] for i in ind_sort]

    return alphas, np.array(mass)


def increment (list):
    list_of_ones=[]
    for item in list:
          list_of_ones.append([1]*len(item))
    incremented_list=[[e1+e2 for e1,e2 in zip(x,y)] for x,y in zip (list,list_of_ones) ] 
    return incremented_list    

def damex(x_norm, R, eps, mu_min):
    x_damex = 1.*(x_norm[np.max(x_norm, axis=1) > R] > R*eps)
    n_extr = np.sum(np.sum(x_damex, axis=1) > 0)
    alphas, mass = damex_0(x_damex)
    n_alphas = np.sum(mass > n_extr * mu_min)
    alphas[:n_alphas]=increment(alphas[:n_alphas])
    return alphas[:n_alphas], mass[:n_alphas]


def list_to_dict_size(list_alphas):
    alphas_dict = {s: [] for s in range(2, max(map(len, list_alphas))+1)}
    for alpha in list_alphas:
        alphas_dict[len(alpha)].append(alpha)

    return alphas_dict