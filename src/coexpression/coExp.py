from typing import Callable
import src.coexpression.pearson as pearson
import src.coexpression.distance as distance
import src.coexpression.rdc as rdc

def pearson_correlation(M: list[list[float]]) -> list[list[float]]:
    '''
    Computes pearson correlation between every pair of rows of a given mxn matrix in O(m^2n)
    '''
    correlation(M, pearson.correlation)

def distance_correlation(M: list[list[float]]) -> list[list[float]]:
    '''
    Computes distance correlation between every pair of rows of a given mxn matrix in O(m^2n^2)
    '''
    correlation(M, distance.correlation)

def rdc_correlation(M: list[list[float]]) -> list[list[float]]:
    '''
    Computes randomized dependence coefficient between every pair of rows of a given mxn matrix in O(m^2nlgn)
    '''
    correlation(M, rdc.correlation)

def correlation(M: list[list[float]], f: Callable[[list[float], list[float]], float]) -> list[list[float]]:
    '''
    Transforms a coexpression matrix M of size mxn into a correlation matrix M' of size mxm.

    m denotes the number of genes, and n denotes the number of samples.
    '''
    m = len(M)

    if m <= 0:
        return [[]]

    C = [[0 for _ in range(m)] for _ in range(m)]

    for i in range(m):
        for j in range(m):
            if j < i:
                C[i][j] = C[j][i]
            elif i == j:
                C[i][j] = 1.0
            else:
                C[i][j] = round(f(M[i], M[j]), 4)

    return C
