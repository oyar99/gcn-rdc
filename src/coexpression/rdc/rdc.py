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

def random_nonlinear_projections():
    pass


def correlation(X: list[float], Y: list[float]) -> float:
    '''
    Computes the randomized dependance coefficient between vectors X and Y of length n in O(nlgn) time
    '''

    return .0

def correlationM(M: list[list[float]]) -> list[list[float]]:
    '''
    Computes the randomized dependence coefficient between all vectors M[i] and M[j] of length n in O(m^2nlgn) time for matrix M of size mxn
    '''
    m = len(M)

    C = [[0 for _ in range(m)] for _ in range(m)]

    # Pre-compute copula transformations for all vectors M[i]
    copulas = [[[] for _ in range(m)] for _ in range(m)]

    for i in range(m):
        copulas[i] = copula(M[i])

    return C
