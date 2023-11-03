from typing import Callable
import src.coexpression.pearson.pearson as pearson
import src.coexpression.distance.distance as distance
import src.coexpression.rdc.rdc as rdc

def pearson_correlation(M: list[list[float]]) -> list[list[float]]:
    '''
    Computes pearson correlation between every pair of rows of a given mxn matrix in O(m^2n)
    '''
    return correlation(M, pearson.correlationM)

def distance_correlation(M: list[list[float]]) -> list[list[float]]:
    '''
    Computes distance correlation between every pair of rows of a given mxn matrix in O(m^2n^2)
    '''
    return correlation(M, distance.correlationM)

def rdc_correlation(M: list[list[float]]) -> list[list[float]]:
    '''
    Computes randomized dependence coefficient between every pair of rows of a given mxn matrix in O(m^2nlgn)
    '''
    return correlation(M, rdc.correlation)

def correlation(M: list[list[float]], f: Callable[[list[list[float]]], float]) -> list[list[float]]:
    '''
    Transforms a coexpression matrix M of size mxn into a correlation matrix M' of size mxm.

    m denotes the number of genes, and n denotes the number of samples.
    '''
    return f(M)
