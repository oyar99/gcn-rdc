class GCN:
    '''
    GCN - Gene Co-expression Network with utility methods for
    statistical analysis.
    '''
    adj_M: list[list[bool]]
    n: int
    m: int

    def __init__(self, M: list[list[float]], t: float):
        '''
        Constructs a GCN from a correlation matrix.
        For all pair of genes, there is an edge between them iff
        their correlation is greater than or equal to t.
        '''
        n = len(M)

        self.n = n
        self.m = 0
        self.adj_M = [[0 for _ in range(n)] for _ in range(n)]

        if n <= 0:
            return
        
        for i in range(n):
            for j in range(len(M[i])):
                self.adj_M[i][j] = abs(M[i][j]) >= t and i != j
                if self.adj_M[i][j]:
                    self.m = self.m + 1

        # Avoid double counting edges
        self.m = self.m / 2

    def degrees(self) -> list[int]:
        degrees = [0 for _ in range(self.n)]
        
        for i in range(len(self.adj_M)):
            degree = 0
            for j in range(len(self.adj_M[i])):
                if self.adj_M[i][j]:
                    degree = degree + 1
            degrees[i] = degree

        return degrees

    def degree_distribution(self) -> list[int]:
        distr = [0 for _ in range(self.n)]

        degrees = self.degrees()

        for d in degrees:
            distr[d] = distr[d] + 1

        return distr
    
    def density(self) -> float:
        return 2 * self.m / (self.n * (self.n - 1))

    def clustering_coefficient(self) -> list[float]:
        clust_coeff = [.0 for _ in range(self.n)]
        degrees = self.degrees()

        for i in range(self.n):
            # Number of adjacent nodes of i that 
            # are connected between other adjacent nodes of i
            neigh_degree = 0
            for j in range(self.n):
                if self.adj_M[i][j]:
                    for k in range(self.n):
                        if self.adj_M[j][k] and self.adj_M[i][k]:
                            neigh_degree = neigh_degree + 1

            # Avoid double counting edges
            neigh_degree = neigh_degree / 2

            if degrees[i] > 1:
                clust_coeff[i] = 2 * neigh_degree / (degrees[i] * (degrees[i] - 1))

        return clust_coeff
        
    
    def spectral_clustering(self) -> list[float]:
        spectral = [.0 for _ in range(self.n)]

        degrees = self.degrees()
        clust_coeff = self.clustering_coefficient()

        for d in range(self.n):
            avg_clust = .0
            k = 0
            for i in range(self.n):
                if degrees[i] == d:
                    avg_clust = avg_clust + clust_coeff[i]
                    k = k + 1

            if k > 0:
                avg_clust = avg_clust / k

            spectral[d] = avg_clust

        return spectral

    def serialize(self, labels: list[str]) -> list[list[str]]:
        M = []

        for i in range(len(self.adj_M)):
            for j in range(i + 1, len(self.adj_M[i])):
                if self.adj_M[i][j]:
                    M.append([labels[i], labels[j]])

        return M