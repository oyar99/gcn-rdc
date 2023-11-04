from typing import Callable
from sklearn.cross_decomposition import CCA
from scipy import stats
import src.coexpression.pearson.pearson as pearson
import numpy as np
import math

def transpose_vector(X: list[float]) -> list[list[float]]:
    n = len(X)

    T = [[0 for _ in range(1)] for _ in range(n)]

    for i in range(n):
        T[i][0] = X[i]

    return T

def vector_to_matrix(X: list[float], r: int) -> list[list[float]]:
    '''
    Given a vector X of length n, returns a matrix of length rxn where each
    row is a copy of X
    '''
    n = len(X)

    Y = [[0 for _ in range(n)] for _ in range(r)]

    for i in range(r):
        for j in range(n):
            Y[i][j] = X[j]

    return Y

def mult_matrix_by_vector(X: list[list[float]], Y: list[float]) -> list[list[float]]:
    '''
    Multiplies a nx1 matrix by a vector of length k
    '''
    n = len(X)
    k = len(Y)

    XY = [[0 for _ in range(k)] for _ in range(n)]

    for i in range(n):
        for j in range(k):
            XY[i][j] = X[i][0] * Y[j]

    return XY

def add_matrices(X: list[list[float]], Y: list[list[float]], n: int, k: int) -> list[list[float]]:
    XY = [[0 for _ in range(k)] for _ in range(n)]

    for i in range(n):
        for j in range(k):
            XY[i][j] = X[i][j] + Y[i][j]

    return XY

def concat_cols(X: list[list[float]], Y: list[list[float]], k: int) -> list[list[float]]:
    n = len(X)

    XY = [[0 for _ in range(2 * k)] for _ in range(n)]

    for i in range(n):
        for j in range(2 * k):
            XY[i][j] = X[i][j] if j < k else Y[i][j - k]

    return XY

def apply_non_linear_function(X: list[list[float]], f: Callable[[float], float] = math.sin) -> list[list[float]]:
    for i in range(len(X)):
        for j in range(len(X[i])):
            X[i][j] = f(X[i][j])

    return X

def copula(X: list[float]) -> list[float]:
    '''
    Computes the empirical copula transformation of a 1-dimensional random variable sample by using
    the cumulative distribution function defined as

    P(x) := 1/n Sum from i=1 to n -> I(Xi <= x)
    '''
    n = len(X)

    x_tuple = [(X[i], i) for i in range(n)]

    x_tuple.sort(key=lambda x: x[0])

    R = [0 for _ in range(n)]
    r = 1

    for (_, i) in x_tuple:
        R[i] = r / n
        r = r + 1

    return R

def random_normal_dist(n: int, s: float) -> list[float]:
    '''
    Generates a random sample of n values from the normal distribution with mean 0 and variance s
    '''
    return np.random.normal(0, s, n).tolist()

def random_uniform_dist(n: int) -> list[float]:
    '''
    Generates a random sample of n values from a uniform distribution in the interval [-PI, PI)
    '''
    return np.random.uniform(-math.pi, math.pi, n).tolist()

def random_nonlinear_projections(A: list[float], k: int, s: float) -> list[list[float]]:
    '''
    Generates random non-linear projections to project a given vector X. The projection is defined by

    P(X, k, s) := 
    
        f(w0x0 + b0) ... f(wkx1 + bk)
        ...
        f(w0xn + b0) ... f(wkxn + bk)

    where f is a non-linear function that will be applied to the random projection.
    Commonly used functions include sine, cosine, parabolas, sinusoids, etc.

    We will use a combination of sine and cosine transformations.

    w denotes a vector of k values drawn from a normal distribution.
    b denotes a vector of k values drawn from a uniform distribution.
    '''
    n = len(A)

    # Let us transpose X so we can then multiply it by vector W
    X = transpose_vector(A)

    W = random_normal_dist(k, s)
    B = vector_to_matrix(random_uniform_dist(k), n)

    # Let us compute XW
    XW_B = add_matrices(mult_matrix_by_vector(X, W), B, n, k)

    return concat_cols(
        apply_non_linear_function(XW_B, math.cos), 
        apply_non_linear_function(XW_B, math.sin), k)

def canonical_correlation_analysis(X: list[list[float]], Y: list[list[float]]) -> float:
    cca = CCA(n_components=1)
    cca.fit(X, Y)

    X_corr, Y_corr = cca.transform(X, Y)
    X_corr = np.hstack(X_corr).tolist()
    Y_corr = np.hstack(Y_corr).tolist()

    X_corr_mean = pearson.mean(X_corr)
    Y_corr_mean = pearson.mean(Y_corr)

    return pearson.correlation(X_corr, Y_corr, pearson.std(X_corr, X_corr_mean), pearson.std(Y_corr, Y_corr_mean), X_corr_mean, Y_corr_mean)

def correlation(X: list[float], Y: list[float], k: int, s: float) -> float:
    '''
    Computes the randomized dependance coefficient between vectors X and Y of length n in O(nlgn) time
    '''

    # Compute random non-linear projections of X and Y in O(nk) time which reduces to O(n) for small k
    fx = random_nonlinear_projections(X, k, s)
    fy = random_nonlinear_projections(Y, k, s)

    # Obtain the maximum correlation coefficients over the projections
    return canonical_correlation_analysis(fx, fy)

def correlationM(M: list[list[float]]) -> list[list[float]]:
    '''
    Computes the randomized dependence coefficient between all vectors M[i] and M[j] of length n in O(m^2nlgn) time for matrix M of size mxn
    '''
    m = len(M)

    C = [[0 for _ in range(m)] for _ in range(m)]

    # Pre-compute copula transformations for all vectors M[i] in O(mnlgn) time
    copulas = [[[] for _ in range(m)] for _ in range(m)]

    for i in range(m):
        copulas[i] = copula(M[i])

    '''
    The coefficient should in practice be very insensitive to the params T, K and S
    '''
    # Let T denote a constant number of times we will compute the randomized dependence coefficient for stability
    T = 10
    # Let K denote the number of random number projections used
    K = 5
    # Let S denote the variance of the normal distribution used for the random number projections
    S = 0.02

    for i in range(m):
        for j in range(m):
            if j < i:
                C[i][j] = C[j][i]
            elif i == j:
                C[i][j] = 1.0
            else:
                corrs = [0 for _ in range(T)]
                for k in range(T):
                    corrs[k] = round(correlation(copulas[i], copulas[j], K, S), 4)

                # Choose the median across all K runs
                corrs.sort()
                C[i][j] = corrs[int(T/2)] if T % 2 == 0 else corrs[int(T/2)] + corrs[int((T-1)/2)]

    return C

'''
Numpy Based Implementation
'''
def correlation_numpy(M: np.ndarray) -> np.ndarray:
    # Pre-compute all copula-transformations for each row
    X = np.apply_along_axis(copulas_numpy, 1, M)

    dcorrelation_inner = lambda Y : np.apply_along_axis(rdc_numpy, 1, X, Y)

    return np.round(np.apply_along_axis(dcorrelation_inner, 1, X), 4)

def rdc_numpy(X: np.ndarray, Y: np.ndarray) -> float:
    '''
    Given the copula transformations of X and Y, project them with non-linear functions
    and use canonical correlation analysis to get the correlation between the original
    vectors X and Y
    '''
    n = len(X)
    # Let T denote a constant number of times we will compute the randomized dependence coefficient for stability
    T = 5
    # Let K denote the number of random number projections used
    K = 4
    # Let S denote the variance of the normal distribution used for the random number projections
    S = 0.15

    A = []

    for _ in range(T):
        x = X[:, np.newaxis]
        Wx = random_normal_dist_numpy(K, S).reshape(1, K)
        Bx = random_uniform_dist_numpy(K).reshape(1, K)

        Xw = np.add(np.dot(x, Wx), np.repeat(Bx, n, axis=0))
        fx = np.concatenate((np.cos(Xw), np.sin(Xw)), axis=1)

        y = Y[:, np.newaxis]
        Wy = random_normal_dist_numpy(K, S).reshape(1, K)
        By = random_uniform_dist_numpy(K).reshape(1, K)

        Yw = np.add(np.dot(y, Wy), np.repeat(By, n, axis=0))
        fy = np.concatenate((np.cos(Yw), np.sin(Yw)), axis=1)

        corr = round(canonical_correlation_analysis_numpy(fx, fy), 4)

        A.append(corr)

    return np.median(np.array([A]))

def copulas_numpy(X: np.ndarray) -> np.ndarray:
    ranks = np.add(X.argsort().argsort(), 1)
    
    return np.multiply(ranks, 1/float(X.shape[0]))

def random_normal_dist_numpy(n: int, s: float) -> np.ndarray:
    '''
    Generates a random sample of n values from the normal distribution with mean 0 and variance s
    '''
    return np.random.normal(0, s, n)

def random_uniform_dist_numpy(n: int) -> np.ndarray:
    '''
    Generates a random sample of n values from a uniform distribution in the interval [-PI, PI)
    '''
    return np.random.uniform(-math.pi, math.pi, n)

def canonical_correlation_analysis_numpy(X: np.ndarray, Y: np.ndarray) -> float:
    cca = CCA(n_components=1)
    cca.fit(X, Y)

    X_corr, Y_corr = cca.transform(X, Y)

    corr, _ = stats.pearsonr(np.hstack(X_corr), np.hstack(Y_corr))

    return corr
