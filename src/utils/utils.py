import pandas as pd

def readCoexpressionFileAsCsv(filepath: str) -> [list[list[float]], list[str]]:
    '''
    Reads a microarray csv file that contains expression data.
    Rows should correspond to genes and columns should have expression data for each gene.

    columns that do not contain expression information are ignored.

    Returns a pandas data frame containing the expression data along with a list of genes in the data
    '''
    drop_cols = {'gene_symbol', 'LocusLinkID', 'ProteomeID', 'cytogeneticLoc', 'CHROMOSOME', 'StartPosition', 'EndPosition'}
    cols = list(pd.read_csv(filepath, nrows=1))

    df = pd.read_csv(filepath, usecols=[i for i in cols if i not in drop_cols])

    if len(df.columns) == 0:
        return [[], []]

    genes = df.iloc[:, 0].to_list()

    np = df.drop(df.columns[0], axis=1).to_numpy()

    M = np.tolist()

    return [M, genes]

def saveFileContent(content: str, filepath: str) -> None:
    with open(filepath, 'w') as f:
        print(content, file=f)

def saveList(L: list[float], filepath: str, labels: list[str] = None) -> None:
    df = pd.DataFrame(L, labels)

    df.to_csv(filepath, sep='\t', header=False)

def saveMatrix(M: list[list[float]], filepath: str, labels: list[str] = None) -> None:
    df = pd.DataFrame(M, labels, labels)

    df.to_csv(filepath, sep='\t', index=labels != None, header=labels != None)


