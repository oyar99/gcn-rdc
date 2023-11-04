import numpy as np

def mean(X: list[float]) -> float:
    sum = .0

    for x in X:
        sum = sum + x

    return sum/len(X)

def covariance(X: list[float], Y: list[float], meanX: float, meanY: float) -> float:
    n = len(X)

    s = .0

    for i in range(n):
        s = s + ((X[i] - meanX) * (Y[i] - meanY))

    return s

def std(X: list[float], meanX: float) -> float:
    s = .0

    for x in X:
        s = s + ((x - meanX) * (x - meanX))

    return s ** .5

def correlation(X: list[float], Y: list[float], stdX: float, stdY: float, meanX: float, meanY: float) -> float:
    '''
    Computes the pearson correlation between vectors X and Y of length n in O(n) time
    '''
    return covariance(X, Y, meanX, meanY)/(stdX * stdY)

def correlationM(M: list[list[float]]) -> list[list[float]]:
    '''
    Computes the pearson correlation between all vectors M[i] and M[j] of length n in O(m^2n) time for matrix M of size mxn
    '''
    m = len(M)

    C = [[0 for _ in range(m)] for _ in range(m)]

    # Pre-compute standard deviations for all vectors M[i]
    S = [[0 for _ in range(m)] for _ in range(m)]
    # Pre-compute means for all vectors M[i]
    A = [[0 for _ in range(m)] for _ in range(m)]

    for i in range(m):
        A[i] = mean(M[i])
        S[i] = std(M[i], A[i])

    for i in range(m):
        for j in range(m):
            if j < i:
                C[i][j] = C[j][i]
            elif i == j:
                C[i][j] = 1.0
            else:
                C[i][j] = round(correlation(M[i], M[j], S[i], S[j], A[i], A[j]), 4)

    return C

'''
Numpy Based Implementation
'''
def correlation_numpy(M: np.ndarray) -> np.ndarray:
    return np.round(np.corrcoef(M), 4)