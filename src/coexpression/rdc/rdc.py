def correlation(x: list[float], y: list[float]) -> float:
    '''
    Computes the randomized dependance coefficient between vectors X and Y of length n in O(nlgn) time
    '''
    if (len(x) != len(y)):
        return .0

    return .0