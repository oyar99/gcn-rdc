import argparse
import numpy as np
import src.gcn.gcn as gcn
import src.coexpression.coExp as coExp
import src.utils.utils as utils

def get_args():
    '''
    Get program arguments specified by the user
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--a', '--analysis', action=argparse.BooleanOptionalAction, help=
                        '''
                        whether the GCN analysis should be output
                        ''')
    parser.add_argument('-t', '--threshold', type=float, default=.75, help=
                        '''
                        correlation threshold. It determines whether two genes are
                        connected or not based on its correlation value
                        ''')
    parser.add_argument('--p', '--performance', action=argparse.BooleanOptionalAction, help=
                        '''
                        whether to use a high performance implementation for correlation measures that heavily use numpy
                        ''')

    requiredNamedArgs = parser.add_argument_group('required named arguments')
    requiredNamedArgs.add_argument('-i', '--input', type=str, required=True, help='Input csv file with coexpression data')
    requiredNamedArgs.add_argument('-cr', '--correlation', type=str, choices=['pearson', 'distance', 'sdistance', 'rdc'],
                                   required=True, help='Correlation measure to use')

    args = parser.parse_args()

    return args

def correlation_unsafe(M: list[list[float]], corr: str) -> list[list[float]]:
    return (coExp.pearson_correlation_unsafe(M) if corr == 'pearson' else
         coExp.distance_correlation_unsafe(M) if corr == 'distance' else
         coExp.signed_distance_correlation_unsafe(M) if corr == 'sdistance' else
         coExp.rdc_correlation_unsafe(M) if corr == 'rdc' else [])

def correlation_numpy(M: list[list[float]], corr: str) -> np.ndarray:
    return (coExp.pearson_correlation(np.array(M)) if corr == 'pearson' else
         coExp.distance_correlation(np.array(M)) if corr == 'distance' else
         coExp.signed_distance_correlation(np.array(M)) if corr == 'sdistance' else
         coExp.rdc_correlation(np.array(M)) if corr == 'rdc' else [])

def main():
    '''
    GCN Entry Point
    '''
    args = get_args()

    # Read coexpression data and transform into matrix
    [M, genes] = utils.read_coexpression_file_as_csv(args.input)
    corr = args.correlation

    C = correlation_unsafe(M, corr) if not args.p else correlation_numpy(M, corr).tolist()

    # Output the pearson correlation matrix with gene labels
    utils.save_matrix(C, 'correlation.csv', labels=genes)

    network = gcn.GCN(C, args.threshold)

    # Serializes and outputs the network to a file where each row corresponds to an edge
    serializedNetwork = network.serialize(genes)
    utils.save_matrix(serializedNetwork, 'network.csv')

    if args.a:
        '''
        Output analysis information
        '''
        degrees = network.degrees()

        utils.save_list(degrees, 'degrees.csv', labels=genes)

        degrees_dist = network.degree_distribution()

        utils.save_list(degrees_dist, 'degree_dist.csv')

        density = network.density()

        utils.save_file_content(f'Network Density: {density}', 'density.txt')

        clust_coeff = network.clustering_coefficient()

        utils.save_list(clust_coeff, 'clustering_coeff.csv', labels=genes)

        spectral = network.spectral_clustering()

        utils.save_list(spectral, 'spectral.csv')

if __name__ == "__main__":
    main()
