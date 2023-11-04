from typing import Callable
import numpy as np
import src.coexpression.pearson.pearson as pearson
import src.coexpression.distance.distance as distance
import src.coexpression.rdc.rdc as rdc

def pearson_correlation_unsafe(M: list[list[float]]) -> list[list[float]]:
    '''
    Computes pearson correlation between every pair of rows of a given mxn matrix in O(m^2n)
    '''
    return correlation_unsafe(M, pearson.correlationM)

def pearson_correlation(np: np.ndarray) -> np.ndarray:
    '''
    Computes pearson correlation between every pair of rows of a given mxn matrix in O(m^2n)
    '''
    return correlation(np, pearson.correlation_numpy)

def distance_correlation_unsafe(M: list[list[float]]) -> list[list[float]]:
    '''
    Computes distance correlation between every pair of rows of a given mxn matrix in O(m^2n^2)
    '''
    return correlation_unsafe(M, distance.correlationM)

def distance_correlation(np: np.ndarray) -> np.ndarray:
    '''
    Computes distance correlation between every pair of rows of a given mxn matrix in O(m^2n^2)
    '''
    return correlation(np, distance.correlation_numpy)

def signed_distance_correlation_unsafe(M: list[list[float]]) -> list[list[float]]:
    '''
    Computes signed distance correlation between every pair of rows of a given mxn matrix in O(m^2n^2)
    '''
    return correlation_unsafe(M, distance.signedCorrelationM)

def signed_distance_correlation(np: np.ndarray) -> np.ndarray:
    '''
    Computes signed distance correlation between every pair of rows of a given mxn matrix in O(m^2n^2)
    '''
    return correlation(np, distance.signed_correlation_numpy)

def rdc_correlation_unsafe(M: list[list[float]]) -> list[list[float]]:
    '''
    Computes randomized dependence coefficient between every pair of rows of a given mxn matrix in O(m^2nlgn)
    '''
    return correlation_unsafe(M, rdc.correlationM)

def rdc_correlation(np: np.ndarray) -> np.ndarray:
    '''
    Computes randomized dependence coefficient between every pair of rows of a given mxn matrix in O(m^2nlgn)
    '''
    return correlation(np, rdc.correlation_numpy)

def correlation_unsafe(M: list[list[float]], f: Callable[[list[list[float]]], float]) -> list[list[float]]:
    '''
    Transforms a coexpression matrix M of size mxn into a correlation matrix M' of size mxm.

    m denotes the number of genes, and n denotes the number of samples.
    '''
    return f(M)

def correlation(np: np.ndarray, f: Callable[[np.ndarray], float]) -> np.ndarray:
    '''
    Transforms a coexpression matrix M of size mxn into a correlation matrix M' of size mxm.

    m denotes the number of genes, and n denotes the number of samples.
    '''
    return f(np)