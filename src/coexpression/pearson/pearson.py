def mean(X: list[float]) -> float:
    sum = .0

    for x in X:
        sum = sum + x

    return sum/len(X)

def covariance(X: list[float], Y: list[float]) -> float:
    n = len(X)

    meanX = mean(X)
    meanY = mean(Y)

    s = .0

    for i in range(n):
        s = s + ((X[i] - meanX) * (Y[i] - meanY))

    return s/(n-1)

def std(X: list[float]) -> float:
    n = len(X)

    meanX = mean(X)

    s = .0

    for x in X:
        s = s + ((x - meanX) * (x - meanX))

    return (s/(n-1)) ** .5

def correlation(X: list[float], Y: list[float]) -> float:
    '''
    Computes the pearson correlation between vectors X and Y of length n in O(n) time
    '''
    if (len(X) != len(Y)):
        return .0

    return covariance(X, Y)/(std(X) * std(Y))