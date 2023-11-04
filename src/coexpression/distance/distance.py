import numpy as np
from dcor import distance_correlation as dcorrelation
import src.coexpression.pearson.pearson as pearson

def pairwise_distances(X: list[float]) -> list[list[float]]:
    n = len(X)

    D = [[0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(n):
            D[i][j] = abs(X[i] - X[j])

    return D

def doubly_centered_distances(X: list[list[float]]) -> list[list[float]]:
    n = len(X)

    row_means = [0 for _ in range(n)]

    for i in range(n):
        mean = .0
        for j in range(n):
            mean = mean + X[i][j]
        row_means[i] = mean/n

    col_means = [0 for _ in range(n)]

    for j in range(n):
        mean = .0
        for i in range(n):
            mean = mean + X[i][j]
        col_means[j] = mean/n
    
    mean = .0

    for i in range(n):
        for j in range(n):
            mean = mean + X[i][j]
    
    mean = mean/(n * n)

    for i in range(n):
        for j in range(n):
            X[i][j] = X[i][j] - row_means[i] - col_means[j] + mean

    return X

def distance_covariance(DX: list[list[float]], DY: list[list[float]]) -> float:
    n = len(DX)

    dcov2 = .0

    for i in range(n):
        for j in range(n):
            dcov2 = dcov2 + (DX[i][j] * DY[i][j])

    return dcov2

def distance_variance(DX: list[float]) -> float:
    n = len(DX)
    s = .0

    for i in range(n):
        for j in range(n):
            s  = s + (DX[i][j] * DX[i][j])

    return s

def correlation(DX: list[list[float]], DY: list[list[float]], dvx: float, dvy: float) -> float:
    '''
    Computes the distance correlation between vectors X and Y of length n in O(n^2) time
    '''
    return (distance_covariance(DX, DY)/((dvx ** .5) * (dvy ** .5))) ** .5

def signedCorrelationM(M: list[list[float]]) -> float:
    '''
    Computes the signed distance correlation between all vectors M[i] and M[j] of length n in O(m^2n^2) time for matrix M of size mxn
    '''
    m = len(M)

    CD = correlationM(M)
    CP = pearson.correlationM(M)

    for i in range(m):
        for j in range(m):
            CD[i][j] = CD[i][j] * (-1 if CP[i][j] < 0 else 1)

    return CD

def correlationM(M: list[list[float]]) -> float:
    '''
    Computes the distance correlation between all vectors M[i] and M[j] of length n in O(m^2n^2) time for matrix M of size mxn
    '''
    m = len(M)

    C = [[0 for _ in range(m)] for _ in range(m)]

    # Pre-compute all doubly centered distances for all vectors M[i]
    D = [[[] for _ in range(m)] for _ in range(m)]
    # Pre-compute all distance variances for all vectors M[i]
    V = [[0 for _ in range(m)] for _ in range(m)]

    for i in range(m):
        D[i] = doubly_centered_distances(pairwise_distances(M[i]))
        V[i] = distance_variance(D[i])

    for i in range(m):
        for j in range(m):
            if j < i:
                C[i][j] = C[j][i]
            elif i == j:
                C[i][j] = 1.0
            else:
                C[i][j] = round(correlation(D[i], D[j], V[i], V[j]), 4)

    return C

'''
Numpy Based Implementation
'''
def correlation_numpy(M: np.ndarray) -> np.ndarray:
    dcorrelation_inner = lambda X : np.apply_along_axis(dcorrelation, 1, M, X)

    return np.round(np.apply_along_axis(dcorrelation_inner, 1, M), 4)

def signed_correlation_numpy(M: np.ndarray) -> np.ndarray:
    C = correlation_numpy(M)
    P = pearson.correlation_numpy(M)

    return np.multiply(C, np.sign(P))