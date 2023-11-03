def pairwise_distance(X: list[float]) -> list[list[float]]:
    n = len(X)

    D = [[0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(n):
            D[i][j] = abs(X[i] - X[j])

    return D

def double_centering(X: list[list[float]]) -> list[list[float]]:
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

def distance(X: list[float], Y: list[float]) -> list[list[float]]:
    DX = pairwise_distance(X)
    DY = pairwise_distance(Y)

    DX = double_centering(DX)
    DY = double_centering(DY)

    




def correlation(x: list[float], y: list[float]) -> float:
    '''
    Computes the distance correlation between vectors X and Y of length n in O(n^2) time
    '''
    if (len(x) != len(y)):
        return .0

    return 0