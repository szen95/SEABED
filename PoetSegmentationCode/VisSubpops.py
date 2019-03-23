
def DrawHierarchicalTree(filename, resultsFolderTag):
    """

            Method to draw hierarchical tree without any outcome shading


    :param filename: string containing the full path and filename that the figure should be written to
    :param numVarsTotal: integer indicating the total number of variables POET digested
    :param resultsFolderTag: the path + expTag to the results folder where the subpop results are retrieved from
    :return: 
    """

    # Imports
    import create_tree_image as cti

    # Get tree parameters
    numVerts, numVarsUsed, leafNames, varsUsed, clinVals, cellLines, gListClosure = RetrieveSubpopResults(resultsFolderTag)
    # junk1, junk2, numVerts, numVarsUsed, leafNames, varsUsed, junk7, junk8, gListInd, gListClosure = RetrieveSubpopResults \
    #     (resultsFolderTag)

    # GEt number of covariates
    numVarsTotal = len(varsUsed[0])

    # ID which subpops have homo vars
    homoVars = list()
    for ii in range(len(numVerts)):
        if numVarsUsed[ii] == numVarsTotal:
            homoVars.append(0)
        else:
            homoVars.append(1)

    # INvoke tree making algorithm
    cti.create_tree(numVerts, leafNames, homoVars, filename)

    return


def DrawHierarchicalTreeWithOutcomes(filename, rfTag1, rfTag2, outcomeStr):
    """

        Method to draw hierarchical tree

    :return: 
    """

    import create_tree_image as cti
    import cPickle

    f = open(rfTag1 + 'CompositeSensitivityVariables', 'r')
    compVar = cPickle.load(f)
    indsSMinMeanCompVar = cPickle.load(f)
    f.close()

    # Retrieve subpop parameters
    numVerts, numVarsUsed, leafNames, varsUsed, clinVals, cellLines, gListClosure = RetrieveSubpopResults(rfTag2)

    avgCompVar = list()
    # for ii in range(compVar[0])

    # GEt number of covariates
    numVarsTotal = len(varsUsed[0])

    # ID which subpops have homo vars
    homoVars = list()
    for ii in range(len(numVerts)):
        if numVarsUsed[ii] == numVarsTotal:
            homoVars.append(0)
        else:
            homoVars.append(1)

    # INvoke tree making algorithm
    # cti.create_tree(numVerts, leafNames, homoVars, filename, outcomeRates, outcomeStr)

    return

def RetrieveSubpopResults(rfTag):
    """

        Method to extract all variables from subpopresults file and return them

    :param rfTag: 
    :return: 
    """

    from cPickle import load


    # Write results to file
    # f = open(localResultsFolder + expTag + curTarg + '_' + 'SubpopResults', 'wb')
    # cPickle.dump(numVerts, f)
    # cPickle.dump(numVarsUsed, f)
    # cPickle.dump(gName, f)
    # cPickle.dump(varsUsed, f)
    # cPickle.dump(clinVals, f)
    # cPickle.dump(cellLines, f)
    # cPickle.dump(gListClosure, f)
    # f.close()



    # Open subpopresults file and load variables
    f = open(rfTag + 'SubpopResults', 'r')
    numVerts = load(f)
    numVarsUsed = load(f)
    gName = load(f)
    varsUsed = load(f)
    clinVals = load(f)
    cellLines = load(f)
    gListClosure = load(f)
    f.close()

    return numVerts, numVarsUsed, gName, varsUsed, clinVals, cellLines, gListClosure


def main():

    import cPickle

    # Define folder in which subpop results exist
    rfTag1 = '/media/nkeshava/HDD1/Output/CellLinePharm/Exp9_Wrapper_Mapk/Results/Exp9_Wrapper_Mapk_124_'
    rfTag2 = '/media/nkeshava/HDD1/Output/CellLinePharm/Exp9_Wrapper_Mapk/Results/124/Exp9_Wrapper_Mapk_124_'

    DrawHierarchicalTreeWithOutcomes('HierarchicalTreeWithOutcomes124.png', rfTag1, rfTag2, 'Comp. Sens.')


    # Define folder in which subpop results exist
    rfTag = '/media/nkeshava/HDD1/Output/CellLinePharm/Exp10_Wrapper_Mapk/Results/125/Exp10_Wrapper_Mapk_125_'

    DrawHierarchicalTree('HierarchicalTree125.png', rfTag)
    x=5

main()

