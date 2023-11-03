def mean(x: list[float]) -> float:
    sum = .0

    for xi in x:
        sum = sum + xi

    return sum/len(x)

def covariance(x: list[float], y: list[float]) -> float:
    if (len(x) != len(y)):
        return .0
    
    n = len(x)

    # compute the mean for each variable
    meanX = mean(x)
    meanY = mean(y)

    s = .0

    for i in range(n):
        s = s + ((x[i] - meanX) * (y[i] - meanY))

    return s/(n-1)

def std(x: list[float]) -> float:
    n = len(x)

    # compute the mean for x
    meanX = mean(x)

    s = .0

    for xi in x:
        s = s + ((xi - meanX) * (xi - meanX))

    return (s/(n-1)) ** .5

def pearsonCorrelation(x: list[float], y: list[float], stdX: float, stdY: float) -> float:
    if (len(x) != len(y)):
        return .0

    return covariance(x, y)/(stdX * stdY)

def correlationM(M: list[list[float]]) -> list[list[float]]:
    '''
    Transforms a coexpression matrix M of size mxn into a correlation matrix M' of size mxm.

    m denotes the number of genes, and n denotes the number of samples.
    '''
    m = len(M)

    if m <= 0:
        return [[]]

    CM = [[0 for _ in range(m)] for _ in range(m)]

    std_T = [0 for _ in range(m)]

    for i in range(m):
        std_T[i] = std(M[i])

    for i in range(m):
        for j in range(m):
            if j < i:
                CM[i][j] = CM[j][i]
            elif i == j:
                CM[i][j] = 1.0
            else:
                CM[i][j] = round(pearsonCorrelation(M[i], M[j], std_T[i], std_T[j]), 4)

    return CM
