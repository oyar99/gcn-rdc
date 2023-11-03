import argparse
import src.gcn.gcn as gcn
import src.coexpression.coExp as coExp
import src.utils.utils as utils

def get_args():
    '''
    Get program arguments specified by the user
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--cr', '--correlation', action=argparse.BooleanOptionalAction, help='whether the correlation matrix should be output')
    parser.add_argument('--a', '--analysis', action=argparse.BooleanOptionalAction, help='whether the GCN analysis should be output')
    parser.add_argument('-t', '--threshold', type=float, default=.75, help=
                        '''
                        correlation threshold. It determines whether two genes are 
                        connected or not based on its correlation value
                        ''')
    
    requiredNamedArgs = parser.add_argument_group('required named arguments')
    requiredNamedArgs.add_argument('-i', '--input', type=str, required=True, help='Input csv file with coexpression data')

    args = parser.parse_args()

    return args

def main():
    '''
    GCN Entry Point
    '''
    args = get_args()

    # Read coexpression data and transform into matrix
    [expressionMatrix, genes] = utils.readCoexpressionFileAsCsv(args.input)

    # For every pair of genes, compute their pearson correlation
    corrMatrix = coExp.correlationM(expressionMatrix)

    if args.cr:
        # Output the pearson correlation matrix with gene labels
        utils.saveMatrix(corrMatrix, 'correlation.csv', labels=genes)

    network = gcn.GCN(corrMatrix, args.threshold)

    # Serializes and outputs the network to a file where each row corresponds to an edge
    serializedNetwork = network.serialize(genes)
    utils.saveMatrix(serializedNetwork, 'network.csv')

    if args.a:
        '''
        Output analysis information
        '''
        degrees = network.degrees()

        utils.saveList(degrees, 'degrees.csv', labels=genes)

        degrees_dist = network.degree_distribution()

        utils.saveList(degrees_dist, 'degree_dist.csv')

        density = network.density()

        utils.saveFileContent(f'Network Density: {density}', 'density.txt')

        clust_coeff = network.clustering_coefficient()

        utils.saveList(clust_coeff, 'clustering_coeff.csv', labels=genes)

        spectral = network.spectral_clustering()

        utils.saveList(spectral, 'spectral.csv')

if __name__ == "__main__":
    main()
