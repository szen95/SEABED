def DefineGlobalCovariateVariables(covs):

    numCovs = len(covs)

    # List names of vars, as in covs
    global allVarsIncl
    allVarsIncl = covs

    # List associated labels, which in this case are the nanes
    global allVarsInclLabels
    allVarsInclLabels = covs

    # If there is only one variable, then the following variables are set automatically
    # If there are more than one in covs, then set based on user parameters in info file

    if numCovs > 1:
        # Indicate type of the covariates
        global allVarsInclType
        allVarsInclType = geneVarsInclType

        global allVarsInclCats
        allVarsInclCats = geneVarsInclCats

    return

def ConvertToBinaryList(x):
    """

        Method to convert an arbitrary length list of integers to a binary valued
        list where the binary value is 1 only if the original value is > 0

    :param x:
    :return:
    """

    binaryList = list()
    for ii in range(len(x)):
        if x[ii] > 0:
            binaryList.append(1)
        else:
            binaryList.append(0)

    return binaryList

def BuildAndSaveMultivariateGraphs(graphName, geneSet='ALL'):
    """

        Method to build a multivariate single-graph, for each individual drug

        This method first aligns the raw query data and formats it for each drug
        of interest using the drug object class; then based on an input value in
        Ntrim, the number of patients for each drug is trimmed to Ntrim; then the
        variables of interest are selected based on the values in inclVars from
        each drug and they are concatenated;  the joint graph is then made

        Outcome data is also extracted for the trimmed population and stored in
        the outcome variable

        This method augments BuildAndSaveJointMultivariateGraphUsingClinVals by adding
        covariates from CDP = conditions, drug exposures, and procedures

    :param Ntrim:
    :return:
    """

    PrintStartTime('Starting to build 1 multivariate drug graph using gene mutations')

    # ONly retain those gene specified by the user; else retain all of them

    if geneSet == 'ALL':
        sys.exit('No ALL option specified')
        # # Define global variables using all genes
        # DefineGlobalCovariateVariables(genes)
        # keepMatchedCellLineMutVals = matchedCellLineMutVals

    else:

        # Iterate thru geneSet and extract appropriate quantities
        geneData = list()
        geneCellLines = list()

        # Iterate thru each gene in the geneset and extract
        for jj in range(len(geneSet)):

            # Check to see whether entry is 'TISSUE'
            if geneSet[jj] != 'TISSUE':

                fileName = expTag + 'public_' + geneVarsSuffix[jj]

                # Read in joint covariates and patient IDs
                f0 = open(savedDataFolder + fileName,'r')
                cellLineMutVal = cPickle.load(f0)
                cellLine = cPickle.load(f0)
                gene = cPickle.load(f0)
                f0.close()

                # Find entry in PivotList(cellLineMutVal)
                geneInd = [ ii for ii,xx in enumerate(gene) if xx == geneSet[jj]]
                if geneInd == []:
                    sys.exit('Could not find ' + geneSet[ii] + ' in ' + fileName)
                else:
                    geneData.append(PivotList(cellLineMutVal)[geneInd[0]])
                    geneCellLines.append(cellLine)

            else:

                # Open up pharmacology data file and get out tissue information
                savedFolder = MakeFolder(opFolder + 'SavedData/')
                f0 = open(savedDataFolder + expTag + 'pharmData','r')
                cellLinePharmData = cPickle.load(f0)
                tissueData = cPickle.load(f0)
                f0.close()

                # Need to turn tissue data into categorical variable
                catKey = list()
                catVal = list()
                catValCtr = 0

                # Iterate thru each tissue and assign categorical integer value
                for kk in range(len(tissueData)):

                    if tissueData[kk] not in catKey:
                        catKey.append(tissueData[kk])
                        catVal.append(catValCtr)
                        catValCtr += 1

                    else:
                        # Find entry in catKey to get index
                        foundInd = [ mm for mm,xx in enumerate(catKey) if xx == tissueData[kk] ]
                        catVal.append(foundInd[0])

                # Append results to overall covariate lists
                geneData.append(catVal)
                geneCellLines.append(cellLinePharmData)

        # Get pharm data cell lines
        savedFolder = MakeFolder(opFolder + 'SavedData/')
        f0 = open(savedDataFolder + expTag + 'pharmData','r')
        cellLinePharmData = cPickle.load(f0)
        f0.close()

        # Get intersection of geneCellLines and cellLinePharmData
        cellLinesX = set(cellLinePharmData)
        for ii in range(len(geneCellLines)):
            cellLinesX.intersection(set(geneCellLines[ii]))

        # Convert back to list and sort
        cellLinesX = sorted(list(cellLinesX))

        # Iterate thru cellLinesX to find out for which cell lines
        # geneData should be extracted

        geneDataX = list()
        for ii in range(len(geneSet)):
            geneDataX.append([])

        for ii in range(len(cellLinesX)):  # current final cell line
            for jj in range(len(geneData)): # current covariate
                curInd = [kk for kk,xx in enumerate(geneCellLines[jj]) if xx == cellLinesX[ii]]
                geneDataX[jj].append(geneData[jj][curInd[0]])

        DefineGlobalCovariateVariables(geneSet)

    # Start making PNG using matched cell line data
    g0, g0EmdMats = MakeMultivariatePatientNetworkGraph(geneDataX, cellLinesX)

    # SAve graph
    # fullGraphName = opFolder + 'SavedData/' + expTag + graphName + '.gt'
    g0.save(graphName, fmt='gt')

    # Save EMD Mats separately for drug0 and drug1
    savedFolder = MakeFolder(opFolder + 'SavedData/')
    f = open(savedFolder + expTag + graphName + '_EmdMats','wb')
    cPickle.dump(g0EmdMats, f)
    f.close()

    PrintEndTime('Finished building and saving 1 multivariate drug graphs using gene mutation values')

    return

def ClassifyCompoundSensitivity(graphPath, drugSensitivityFile, drugNum, suffix, fileStr='Univar'):


    # Save final lists of data for current processing
    f1 = open(graphPath + 'CELLGENEDATA_DRUG'+str(drugNum),'r')
    geneKeepFinal = cPickle.load(f1)
    cellLineKeepFinal = cPickle.load(f1)
    sensValsFinal = cPickle.load(f1)
    binSensValsFinal = cPickle.load(f1)
    f1.close()

    Nsens = sum(binSensValsFinal)
    Nres = len(binSensValsFinal) - Nsens

    # Get list of files in graphPath that end in '.gt'
    # print os.listdir(graphPath)
    graphFileList = list()
    filelist = os.listdir(graphPath)
    for ii in range(len(filelist)):
        if filelist[ii].endswith('.gt'):
            graphFileList.append(filelist[ii])

    graphFileList = sorted(graphFileList)

    # graphFileList = [ xx for ii,xx in enumerate(os.listdir(graphPath)) if any(xx.endswith('.gt')) ]

    # Iterate thru each .gt file and get list of genes for which there are graphs
    geneList = list()
    precision = list()
    recall = list()
    fbeta = list()
    prNK = list()
    reNK = list()
    fbNK = list()
    for ii in range(len(graphFileList)):
        # print ii
        curFile = graphFileList[ii]
        # Find position of '_' in curFile
        # head0, tail0 = os.path.split(curFile)
        # print tail0
        uScoreInd = [jj for jj,xx in enumerate(curFile) if xx == '_']
        # print uScoreInd
        NScore = len(uScoreInd)

        # Extract gene name
        curGene = curFile[0:uScoreInd[NScore-2]]
        geneList.append(curGene)

        # print head0, tail0

        # Open up svd file
        svdName = curFile[0:uScoreInd[NScore-1]] + '_SVD_' + 'DRUG' + str(drugNum)
        f1 = open(graphPath + svdName)
        curLap = cPickle.load(f1)
        curU = cPickle.load(f1)
        # curS = cPickle.load(f1)
        # curV = cPickle.load(f1)
        f1.close()

        # Get fiedler eigenvector for second smallest singular value

        N = len(curU)
        fVec = curU[:,N-2]

        # Convert fVec into binary vector where vals < 0 --> 0 and vals > 0 --> 1
        binFVec = list()
        for jj in range(len(fVec)):
            if fVec[jj] > 0.0:
                binFVec.append(1)
            else:
                binFVec.append(0)

        # Classify vertices in curU based on values in binSensValsFinal
        # Use sklearn.metrics.precision_recall_fscore_support
        # precision, recall, fbeta_score, support = prfs(binSensValsFinal, binFVec,1.0,average='weighted')
        # precision, recall, fbeta_score, support = prfs(binSensValsFinal, binFVec,1.0,average=None)
        sVecs = list()
        for jj in range(len(binSensValsFinal)):
            sVecs.append(binSensValsFinal[jj] + binFVec[jj])
        dVecs = list()
        for jj in range(len(binSensValsFinal)):
            dVecs.append(binSensValsFinal[jj] - binFVec[jj])

        tps = [ jj for jj,xx in enumerate(sVecs) if xx == 2]
        fps = [ jj for jj,xx in enumerate(dVecs) if xx == -1 ]
        fns = [ jj for jj,xx in enumerate(dVecs) if xx == 1 ]
        pr = np.divide( float(len(tps)), float(len(tps) + len(fps)))
        re = np.divide( float(len(tps)), float(len(tps) + len(fns)))
        if pr > 0 and re > 0:
            fb = np.divide(2.0 * pr * re,pr+re)
        else:
            fb = np.float64(0.0)

        curPrecision, curRecall, curFbeta_score, curSupport = prfs(np.array(binSensValsFinal), np.array(binFVec),1.0,pos_label=None,average='weighted')
        precision.append(curPrecision)
        recall.append(curRecall)
        fbeta.append(curFbeta_score)

        prNK.append(pr)
        reNK.append(re)
        fbNK.append(fb)

        # print out
        # print ii, curGene, curPrecision, curRecall, curFbeta_score, curSupport
        # print ii, curGene, pr, re, fb


    # Get gene index for max value of each
    maxPreInd = [ii for ii,xx in enumerate(prNK) if all(xx >= prNK)]
    maxReInd = [ii for ii,xx in enumerate(reNK) if all(xx >= reNK) ]
    maxFInd = [ii for ii,xx in enumerate(fbNK) if all(xx >= fbNK) ]

    # Print out top 10 of precision, recall, and f1
    sPrecision = np.sort(prNK)
    sPrecisionInds = np.argsort(prNK)
    sRecall = np.sort(reNK)
    sRecallInds = np.argsort(reNK)
    sFbeta = np.sort(fbNK)
    sFbetaInds = np.argsort(fbNK)

    print 'TOP PRECISION GENES'
    print '-------------------'
    for ii in range(len(prNK)-1,len(prNK)-11,-1):
        print len(prNK)-ii, sPrecision[ii], geneList[sPrecisionInds[ii]]
    print
    print 'TOP RECALL GENES'
    print '----------------'
    for ii in range(len(reNK)-1,len(reNK)-11,-1):
        print len(reNK)-ii, sRecall[ii], geneList[sRecallInds[ii]]
    print
    print 'TOP FBETA GENES'
    print '----------------'
    for ii in range(len(fbNK)-1,len(fbNK)-11,-1):
        print len(fbNK)-ii, sFbeta[ii], geneList[sFbetaInds[ii]]
    print

    #########################################################
    # FIGURES
    #########################################################

    if 'univar' in lower(fileStr):
        # Make title string
        tStr1 = 'Precision, Recall, and F-Beta for Drug ' + str(drugNum)
        tStr2 = str(Nsens) + ' Sensitive Cell Lines (AUC TH = ' + str(sensAucTh) + '),  ' + str(Nres) + ' Resistant'
        tStr = tStr1 + '\n' + tStr2

        # Plot results
        fig = plt.figure()
        fig.set_size_inches(10,8)

        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)

        # Plot fractional cdf of counts to show counts vs. # of edges
        ax1.plot(range(0,len(prNK)), prNK)
        # Add red asterisks for max values
        for ii in range(len(maxPreInd)):
            ax1.plot(maxPreInd[ii], prNK[maxPreInd[ii]], 'r*')
        ax1.set_ylim([0,1.0])
        ax1.grid()
        ax1.set_xlabel('')
        ax1.set_ylabel('Precision')
        ax1.set_title(tStr)
        ax1.set_xticks(list(range(len(prNK))))
        ax1.set_xticklabels([]*len(prNK))
        ax1.autoscale(axis='x',tight=True)

        ax2.plot(range(0,len(reNK)), reNK)
        # Add red asterisks for max values
        for ii in range(len(maxPreInd)):
            ax2.plot(maxReInd[ii], reNK[maxReInd[ii]], 'r*')
        ax2.set_ylim([0,1.0])
        ax2.grid()
        ax2.set_xlabel('')
        ax2.set_ylabel('Recall')
        ax2.set_title('')
        ax2.set_xticks(list(range(len(recall))))
        ax2.set_xticklabels([]*len(recall))
        ax2.autoscale(axis='x',tight=True)

        ax3.plot(range(0,len(fbNK)), fbNK)
        # Add red asterisks for max values
        for ii in range(len(maxFInd)):
            ax3.plot(maxFInd[ii], fbNK[maxFInd[ii]], 'r*')
        ax3.set_ylim([0,1.0])
        ax3.grid()
        ax3.set_xlabel('')
        ax3.set_ylabel('F-Beta')
        ax3.set_title('')
        ax3.set_xticks(list(range(len(geneList))))
        ax3.set_xticklabels(geneKeepFinal, rotation='vertical')
        ax3.autoscale(axis='x',tight=True)

        plt.savefig(figFolder + fileStr + 'Precision_' + suffix + '.png')
        plt.close()

    # elif 'bivar' in lower(fileStr):
    #
    #     # Create and then populate array with bivar F1 values
    #     recallMatrix = np.zeros((len(geneKeepFinal),len(geneKeepFinal)))
    #
    #     for ii in range(len(geneKeepFinal)):
    #
    #         for ii in range(ii+1,len(geneKeepFinal)):
    #             recallMatrix.itemset((ii,jj), 3)



    return precision, recall, fbeta

def CorrelateBinaryMutationsAcrossCellLines(univarGraphsFolder, mutationDataFile, suffix, drugNum):

    # Open mutation data for drug num
    f0 = open(mutationDataFile, 'r')
    mutVal = cPickle.load(f0)    # cell lines x genes
    cellLine = cPickle.load(f0)
    gene = cPickle.load(f0)
    f0.close()

    # Open up cellgenedata file
    f1 = open(univarGraphsFolder + 'CELLGENEDATA_DRUG'+str(drugNum), 'r')
    geneKeepFinal = sorted(cPickle.load(f1))
    cellLineFinal = cPickle.load(f1)
    sensVals = cPickle.load(f1)
    binSensVals = cPickle.load(f1)
    f1.close()

    # Find mutvalues for the cell lines in cellLineKeep and genes in geneKeepFinal
    cellLineInds = list()
    for ii in range(len(cellLineFinal)):
        cellLineInds.append(cellLine.index(cellLineFinal[ii]))

    geneInds = list()
    for ii in range(len(geneKeepFinal)):
        geneInds.append(gene.index(geneKeepFinal[ii]))

    # Retain elements of mutVal pertaining to cellLineInds, geneInds
    mutValTrunc = list()
    numMutsByCellLine = list()
    for ii in cellLineInds:
        curGeneVals = list()
        for jj in geneInds:
            curGeneVals.append(abs(mutVal[ii][jj]))
        mutValTrunc.append(curGeneVals)
        numMutsByCellLine.append(sum(curGeneVals))

    # Look at problem ax gene x cell Line problem
    mutValGene = PivotList(mutValTrunc)
    numMutsByGene = list()
    for ii in range(len(mutValGene)):
        numMutsByGene.append(sum(mutValGene[ii]))

    # Calculate dot product
    geneCorr = np.zeros((len(geneInds),len(geneInds)))
    for ii in range(len(geneInds)):
        for jj in range(ii+1,len(geneInds)):
            geneCorr.itemset((ii,jj), np.dot(mutValGene[ii],mutValGene[jj]))
            geneCorr.itemset((jj,ii), np.dot(mutValGene[ii],mutValGene[jj]))

    # Plot results
    fig = plt.figure()
    fig.set_size_inches(10,8)
    ax1 = fig.add_subplot(111)

    # Plot fractional cdf of counts to show counts vs. # of edges
    ax1.bar(range(len(numMutsByGene)), numMutsByGene)
    ax1.grid()
    ax1.set_xlabel('Gene')
    ax1.set_ylabel('Counts')
    ax1.set_title('Number of Cell Lines Where ' + suffix + ' Is Present')
    ax1.set_xticks(list(range(len(geneKeepFinal))))
    ax1.set_xticklabels(geneKeepFinal, rotation='vertical')
    ax1.autoscale(axis='x',tight=True)
    plt.savefig(figFolder + 'CellLineCountsByGene_' + suffix + '.png')
    plt.close()

    fig = plt.figure()
    fig.set_size_inches(10,8)
    ax1 = fig.add_subplot(111)
    im1 = ax1.pcolor(geneCorr)
    ax1.set_xticks(list(range(len(geneKeepFinal))))
    ax1.set_xticklabels(geneKeepFinal, rotation='vertical')
    ax1.set_yticks(list(range(len(geneKeepFinal))))
    ax1.set_yticklabels(geneKeepFinal)
    for item in ax1.get_xticklabels()+ax1.get_yticklabels():
        item.set_fontsize(10)
    ax1.autoscale(axis='both',tight=True)
    plt.colorbar(im1)
    plt.savefig(figFolder + 'GeneCorrelationMatrix_' + suffix + '.png')
    plt.close()

    return

def PerformGraphSvd(graphName=opFolder+'SavedData/'+expTag+'Drug0MultivariateGraph.gt'):

    PrintStartTime('Beginning segmentation of ' + graphName)

    # Load graph
    g = load_graph(graphName, fmt='gt')

    # Get weighted adjacency matrix
    # PrintStartTime('Started calculating graph laplacian')
    adj = adjacency(g, g.edge_properties['simProp']).todense()
    lap = PerformLaplacian(adj)


    # lap = laplacian(g, deg='total', weight=g.edge_properties['simProp']).todense()

    # Perform svd of lap
    PrintStartTime('Started SVD of laplacian')
    u,s,v = np.linalg.svd(lap)

    # Save SVD components to file

    # Read in joint covariates and patient IDs
    PrintStartTime('Saving results of segmentation')
    # savedFolder = MakeFolder(opFolder + 'SavedData/')
    f0 = open(graphName[0:-3] + '_SvdComponents','wb')
    # f0 = open(savedFolder + expTag + 'SvdComponents','wb')
    cPickle.dump(lap, f0)
    cPickle.dump(u, f0)
    cPickle.dump(s, f0)
    cPickle.dump(v, f0)
    f0.close()

    return

def CalcBpccUsingSvdForOneGene(graphPath):
    """
        Method to analyze the SVD of the laplacian, which has already been computed
        in a previous step.
    :return:
    """

    # Read in joint covariates and patient IDs
    PrintStartTime('Saving results of segmentation')
    # savedFolder = MakeFolder(opFolder + 'SavedData/')
    f0 = open(graphPath[0:-3] + '_SvdComponents','r')
    cPickle.load(f0)
    u = cPickle.load(f0)
    s = cPickle.load(f0)
    v = cPickle.load(f0)
    f0.close()

    # Load graph
    g = load_graph(graphPath)
    cellLineProp  = g.vertex_properties['idProp']

    # Find member in eah partition, using u1
    N = len(u)
    u1 = u[:,N-2]
    class1 = list()
    class2 = list()

    for ii in range(len(u1)):
        if u1[ii] < 0:
            class1.append(ii)
        else:
            class2.append(ii)

    # Read in pharmacology data to see whether the partitioning enriches
    # the sensitivity wrt/ certain drugs
    #
    # Read in joint covariates and patient IDs
    f0 = open(savedDataFolder + expTag + inputSensitivityFile[0][1],'r')
    cellLineUn = cPickle.load(f0)
    tissueUn = cPickle.load(f0)
    drug1SensUn = cPickle.load(f0)
    drug2SensUn = cPickle.load(f0)
    drug3SensUn = cPickle.load(f0)
    sensDataUn = cPickle.load(f0)
    headerUn = cPickle.load(f0)
    f0.close()

    # Pharmacology data doesn't list the cell lines in the same order as the graph
    # vertices:  hence, reorder the entries in the pharmacology data to match the
    # ordering in cellLineProp

    gCellLine = list()
    for ii in range(g.num_vertices()):
        gCellLine.append(cellLineProp[g.vertex(ii)])

    g2PhInds = list()
    for ii in range(len(gCellLine)):
        g2PhInds.append(cellLineUn.index(gCellLine[ii]))

    # Reorder pharmacology data using ordering in g2PhInds
    cellLine = list()
    tissue = list()
    drug1Sens = list()
    drug2Sens = list()
    drug3Sens = list()
    sensData = list()

    for ii in g2PhInds:
        cellLine.append(cellLineUn[ii])
        tissue.append(tissueUn[ii])
        drug1Sens.append(drug1SensUn[ii])
        drug2Sens.append(drug2SensUn[ii])
        drug3Sens.append(drug3SensUn[ii])
        sensData.append(sensDataUn[ii])

    # Separate drugs by whether their sensitivity (AUC) is > or < 0.80
    lowSens1 = list()
    highSens1 = list()
    na1 = list()
    for ii in range(len(drug1Sens)):
        if drug1Sens[ii] == 'NA':
            na1.append(ii)
        else:
            if float(drug1Sens[ii]) < 0.80:
                highSens1.append(ii)
            else:
                lowSens1.append(ii)

    lowSens2 = list()
    highSens2 = list()
    na2 = list()
    for ii in range(len(drug2Sens)):
        if drug2Sens[ii] == 'NA':
            na2.append(ii)
        else:
            if float(drug2Sens[ii]) < 0.80:
                highSens2.append(ii)
            else:
                lowSens2.append(ii)

    lowSens3 = list()
    highSens3 = list()
    na3 = list()
    for ii in range(len(drug3Sens)):
        if drug3Sens[ii] == 'NA':
            na3.append(ii)
        else:
            if float(drug3Sens[ii]) < 0.80:
                highSens3.append(ii)
            else:
                lowSens3.append(ii)

    # Score performance by finding % of class1 cell lines in lowSens and class2
    # cell lines in highSens

    score1 = list()
    score2 = list()
    score3 = list()
    for ii in range(len(cellLine)):
        # Score drug 1
        if drug1Sens[ii] != 'NA':
            if ii in lowSens1:
                if ii in class1:
                    score1.append(1)
                else:
                    score1.append(-1)
            else:
                if ii in class2:
                    score1.append(2)
                else:
                    score1.append(-2)
        else:
            score1.append(0)

        # Score drug 2
        if drug2Sens[ii] != 'NA':
            if ii in lowSens2:
                if ii in class1:
                    score2.append(1)
                else:
                    score2.append(-1)
            else:
                if ii in class2:
                    score2.append(2)
                else:
                    score2.append(-2)
        else:
            score2.append(0)

        # Score drug 3
        if drug3Sens[ii] != 'NA':
            if ii in lowSens3:
                if ii in class1:
                    score3.append(1)
                else:
                    score3.append(-1)
            else:
                if ii in class2:
                    score3.append(2)
                else:
                    score3.append(-2)
        else:
            score3.append(0)

    # Get indices for right/wrong cell lines for each drug
    drug1Correct1 = len( [ ii for ii,xx in enumerate(score1) if xx == 1 ])
    drug1Wrong1   = len( [ ii for ii,xx in enumerate(score1) if xx == -1 ])
    drug1Correct2 = len( [ ii for ii,xx in enumerate(score1) if xx == 2 ])
    drug1Wrong2   = len( [ ii for ii,xx in enumerate(score1) if xx == -2 ])

    drug2Correct1 = len( [ ii for ii,xx in enumerate(score2) if xx == 1 ])
    drug2Wrong1   = len( [ ii for ii,xx in enumerate(score2) if xx == -1 ])
    drug2Correct2 = len( [ ii for ii,xx in enumerate(score2) if xx == 2 ])
    drug2Wrong2   = len( [ ii for ii,xx in enumerate(score2) if xx == -2 ])

    drug3Correct1 = len( [ ii for ii,xx in enumerate(score3) if xx == 1 ])
    drug3Wrong1   = len( [ ii for ii,xx in enumerate(score3) if xx == -1 ])
    drug3Correct2 = len( [ ii for ii,xx in enumerate(score3) if xx == 2 ])
    drug3Wrong2   = len( [ ii for ii,xx in enumerate(score3) if xx == -2 ])

    drugResults = list()
    drugResults.append([[drug1Correct1, drug1Correct1+drug1Wrong1],[drug1Correct2, drug1Correct2+drug1Wrong2]])
    drugResults.append([[drug2Correct1, drug2Correct1+drug2Wrong1],[drug2Correct2, drug2Correct2+drug2Wrong2]])
    drugResults.append([[drug3Correct1, drug3Correct1+drug3Wrong1],[drug3Correct2, drug3Correct2+drug3Wrong2]])

    # Get balanced PCC
    drugBpcc = list()

    # Do first drug
    bpcc1 = 100 * np.divide(float(drug1Correct1),float(drug1Correct1+drug1Wrong1))
    bpcc2 = 100 * np.divide(float(drug1Correct2),float(drug1Correct2+drug1Wrong2))
    if 0.5*(bpcc1 + bpcc2) >= 0.50:
        drugBpcc.append( 0.50*(bpcc1 + bpcc2))
    else:
        drugBpcc.append( 1.0 - 0.50*(bpcc1 + bpcc2))

    bpcc1 = 100 * np.divide(float(drug2Correct1),float(drug2Correct1+drug2Wrong1))
    bpcc2 = 100 * np.divide(float(drug2Correct2),float(drug2Correct2+drug2Wrong2))
    if 0.5*(bpcc1 + bpcc2) >= 0.50:
        drugBpcc.append( 0.50*(bpcc1 + bpcc2))
    else:
        drugBpcc.append( 1.0 - 0.50*(bpcc1 + bpcc2))

    bpcc1 = 100 * np.divide(float(drug3Correct1),float(drug3Correct1+drug3Wrong1))
    bpcc2 = 100 * np.divide(float(drug3Correct2),float(drug3Correct2+drug3Wrong2))
    if 0.5*(bpcc1 + bpcc2) >= 0.50:
        drugBpcc.append( 0.50*(bpcc1 + bpcc2))
    else:
        drugBpcc.append( 1.0 - 0.50*(bpcc1 + bpcc2))

    # Print out BPCC's
    print
    print 'Results for using gene set = ' + str(geneVarsInclLabels)
    print '-------------------------------------'
    for ii in range(len(drugName)):
        curPcc = list()
        for jj in range(0,2):
            curData = drugResults[ii][jj]
            curPcc.append(100.0 * np.divide(float(curData[0]),float(curData[1])))

        print str(drugResults[ii]), str(curPcc)
        print 'BPCC for ' + drugName[ii] + ' is ' + str(drugBpcc[ii])

    # Visualize cell line
    VisualizeCellLineNetwork(g, highSens3, lowSens3, score3)

    # Analyze subspace power of highSens and lowSens relVecs
    highSensFracPower = GetSubspaceRelationship(g, highSens3, lowSens3)

    return drugBpcc, highSensFracPower

def CalcBpccUsingSvd(geneName, graphFolder, inputSensPh, drugNum, suffix):
    """
        Method to analyze the SVD of the laplacian, which has already been computed
        in a previous step.
    :return:
    """

    # Construct svd file path
    svdPath = graphFolder + expTag + geneName + '_' + suffix + '_' + 'SVD' + '_' + 'DRUG' + str(drugNum)
    # svdPath = graphFolder + geneName + '_' + suffix + '_' + 'SVD' + '_' + 'DRUG' + str(drugNum)

    # Open up cellgenedata file to get drug sensitivity
    f1 = open(graphFolder + '/CELLGENEDATA_DRUG' + str(drugNum),'r')
    geneKeepFinal = cPickle.load(f1)
    cellLineKeepFinal = cPickle.load(f1)
    sensValsFinal = cPickle.load(f1)
    binSensValsFinal = cPickle.load(f1)
    f1.close()

    # Open up svd file
    f1 = open(svdPath,'r')
    curLap = cPickle.load(f1)
    curU = cPickle.load(f1)
    # curS = cPickle.load(f1)
    # curV = cPickle.load(f1)
    f1.close()

    # Get fiedler eigenvector for second smallest singular value
    N = len(curU)
    fVec = curU[:,N-2]

    # Convert fVec into binary vector where vals < 0 --> 0 and vals > 0 --> 1
    binFVec = list()
    for jj in range(len(fVec)):
        if fVec[jj] > 0.0:
            binFVec.append(1)
        else:
            binFVec.append(0)

    # Binary classify the GBM segmentation results against the actual sensitivity results
    class0C, class0W, class1C, class1W = ScoreBinaryClassifier(binFVec,binSensValsFinal)

    # Get balanced PCC
    bpcc0 = np.divide(float(len(class0C)),float(len(class0C)+len(class0W)))
    bpcc1 = np.divide(float(len(class1C)),float(len(class1C)+len(class1W)))
    if 0.5*(bpcc0 + bpcc1) >= 0.50:
        bpcc = 0.50*(bpcc0 + bpcc1)
    else:
        bpcc = 1.0 - 0.50*(bpcc0 + bpcc1)

    return bpcc

def VisualizeBinaryClassification(geneName, graphFolder, inputSensPh, drugNum, suffix):

    # Construct svd file path
    svdPath = graphFolder + expTag + geneName + '_' + suffix + '_' + 'SVD' + '_' + 'DRUG' + str(drugNum)

    # Open up cellgenedata file to get drug sensitivity
    f1 = open(graphFolder + '/CELLGENEDATA_DRUG' + str(drugNum),'r')
    geneKeepFinal = cPickle.load(f1)
    cellLineKeepFinal = cPickle.load(f1)
    sensValsFinal = cPickle.load(f1)
    binSensValsFinal = cPickle.load(f1)
    f1.close()

    # Open up svd file
    f1 = open(svdPath,'r')
    curLap = cPickle.load(f1)
    curU = cPickle.load(f1)
    # curS = cPickle.load(f1)
    # curV = cPickle.load(f1)
    f1.close()

    # Get fiedler eigenvector for second smallest singular value
    N = len(curU)
    fVec = curU[:,N-2]

    # Convert fVec into binary vector where vals < 0 --> 0 and vals > 0 --> 1
    binFVec = list()
    for jj in range(len(fVec)):
        if fVec[jj] > 0.0:
            binFVec.append(1)
        else:
            binFVec.append(0)

    # Binary classify the GBM segmentation results against the actual sensitivity results
    class0Correct, class0Wrong, class1Correct, class1Wrong = ScoreBinaryClassifier(binFVec,binSensValsFinal)

    # Visualize graph of cell lines with coloring based on classification results
    graphName = svdPath.replace('_SVD','') + '.gt'
    g = load_graph(graphName)

    # Create vertex property for truth value based on binSensValsFinal
    truthProp = g.new_vertex_property('bool')
    for ii in range(g.num_vertices()):
        if binSensValsFinal[ii]:
            truthProp[g.vertex(ii)] = 1
        else:
            truthProp[g.vertex(ii)] = 0

    correctProp = g.new_vertex_property('bool')
    for ii in range(g.num_vertices()):
        if ii in class0Correct or ii in class1Correct:
            correctProp[g.vertex(ii)] = 1
        else:
            correctProp[g.vertex(ii)] = 0

    # Get precision information


    precision, recall, fbeta = CalcPrecisionRecallF1NK(binSensValsFinal, binFVec)
    # precision, recall, fbeta, support = prfs(binSensValsFinal, binFVec,1.0,pos_label=None,average='weighted')

    # Create edge filter to eliminate edges and turn off edges except for when ccProp3
    # is 1 for both vertices of the edge
    edgeFilterProp = g.new_edge_property('bool')
    for ii in range(g.num_vertices()):
        for jj in range(g.num_vertices()):
            if ii in class0Wrong and jj in class0Wrong and ii != jj:
                edgeFilterProp[g.edge(g.vertex(ii),g.vertex(jj))]=1
            elif ii in class1Wrong and jj in class1Wrong and ii != jj:
                edgeFilterProp[g.edge(g.vertex(ii),g.vertex(jj))]=1
            else:
                edgeFilterProp[g.edge(g.vertex(ii),g.vertex(jj))]=0

    # Set edge filter of g to values in edgeFilterProp
    g.set_edge_filter(edgeFilterProp)

    cellLineProp = CreatePatientDataVertexProperty(g,cellLineKeepFinal,'string')

    # Make network image using
    graphName = figFolder + geneName + '_' + suffix + '_' + 'GBMSEG_' + 'DRUG' + str(drugNum) + '.png'
    graph_draw(g, vertex_text_color='r', vertex_fill_color=truthProp, vertex_size=10, \
               vertex_text = cellLineProp, vertex_font_size=12, vertex_halo=True, vertex_halo_color = correctProp, output_size=(4000,4000), output=graphName)
    g.set_edge_filter(None)

    return [ class0Correct, class0Wrong, class1Correct, class1Wrong ], [precision, recall, fbeta ]

def ScoreBinaryClassifier(xPred, xTrue):
    """

        Method to generate sample-level classifier accuracy results from two
        binary-valued vectors

    :param xPred:
    :param xTrue:
    :return:
    """

    class0Correct = list()
    class0Wrong = list()
    class1Correct = list()
    class1Wrong = list()

    # Iterate thru each value in xTrue and determine if it is classified correctly
    for ii in range(len(xTrue)):
        if xTrue[ii] == 1:
            if xPred[ii] == 1:
                class1Correct.append(ii)
            else:
                class1Wrong.append(ii)
        elif xTrue[ii] == 0:
            if xPred[ii] == 1:
                class0Wrong.append(ii)
            else:
                class0Correct.append(ii)

    return class0Correct, class0Wrong, class1Correct, class1Wrong

def CalcPrecisionRecallF1NK(binTruth, binPred):
    """
        Method to calclulate precision, recall, f1

    :param binSensValsFinal: binary list of
    :param binFVec:
    :return:
    """

    # Make list that is sum of inputs
    sVecs = list()
    for jj in range(len(binTruth)):
        sVecs.append(binTruth[jj] + binPred[jj])

    # Make list that is difference of inputs
    dVecs = list()
    for jj in range(len(binTruth)):
        dVecs.append(binTruth[jj] - binPred[jj])

    tps = [ jj for jj,xx in enumerate(sVecs) if xx == 2]
    fps = [ jj for jj,xx in enumerate(dVecs) if xx == -1 ]
    fns = [ jj for jj,xx in enumerate(dVecs) if xx == 1 ]
    pr = np.divide( float(len(tps)), float(len(tps) + len(fps)))
    re = np.divide( float(len(tps)), float(len(tps) + len(fns)))
    if pr > 0 and re > 0:
        fb = np.divide(2.0 * pr * re,pr+re)
    else:
        fb = np.float64(0.0)


    return pr, re, fb

def FindCellLineDifferences(cellLineActiveAlts, binPart):

    numAlts = len(cellLineActiveAlts[0])
    num0 = [ ii for ii,xx in enumerate(binPart) if xx == 0 ]
    num1 = [ ii for ii,xx in enumerate(binPart) if xx == 1 ]

    # Iterate thru entire list of different types of alts in each cell line and
    # collect the total list of cell lines available in each
    cellLineMaster = list()
    for ii in range(numAlts):
        cellLineMaster.append([])

    # Iterate thru each type of alteration
    altParts = list()
    for ii in range(numAlts):
        # Iterate thru each cell line and add gene list to one of 2 lists,
        # based on the value in binPart
        list0 = list()
        list1 = list()
        for jj in range(len(binPart)):
            if binPart[jj] == 0:
                list0 = list0 + cellLineActiveAlts[jj][ii]
            elif binPart[jj] == 1:
                list1 = list1 + cellLineActiveAlts[jj][ii]
        altParts.append([list0, list1])

    # Perform set intersections of genes in altParts to identify
    # genes uniquely in one partition, but not the other
    altGeneDiffs = list()
    for ii in range(numAlts):
        curGeneDiff = list()
        # Convert gene lists for current genetic alteration into lists
        sAltParts0 = set(altParts[ii][0])
        sAltParts1 = set(altParts[ii][1])

        # Get differences in both directions to find out what is in each
        # partition class that is not in the other
        curGeneDiff.append(sAltParts0.difference(sAltParts1))
        curGeneDiff.append(sAltParts1.difference(sAltParts0))

        altGeneDiffs.append(curGeneDiff)

    return altGeneDiffs, altParts

def VisualizeCellLineNetwork(g, highSens3, lowSens3, score3):


    cellLineProp  = g.vertex_properties['idProp']


    # Visualize graph network of cell lines colored by sensitivity
    # and whether they were classified correctly or not
    drugSensProp3 = g.new_vertex_property('float')
    ccProp3 = g.new_vertex_property('float')

    # Iterate thru each vertex
    for ii in range(g.num_vertices()):
        # Set drugSensProp3 to 1 if ii in highSens3
        if ii in highSens3:
            drugSensProp3[g.vertex(ii)] = 1.0
        else:
            drugSensProp3[g.vertex(ii)] = 0.0
        # Set ccProp3 to 1 if score > 0
        if score3[ii] > 0:
            ccProp3[g.vertex(ii)] = 1
        else:
            ccProp3[g.vertex(ii)] = 0

    # Create edge filter to eliminate edges and turn off edges except for when ccProp3
    # is 1 for both vertices of the edge
    edgeFilterProp = g.new_edge_property('bool')
    for ii in range(g.num_vertices()):
        for jj in range(g.num_vertices()):
            if drugSensProp3[g.vertex(ii)] == 15 and drugSensProp3[g.vertex(jj)] == 1:
                edgeFilterProp[g.edge(g.vertex(ii),g.vertex(jj))]=1
            else:
                edgeFilterProp[g.edge(g.vertex(ii),g.vertex(jj))]=0

    # Set edge filter of g to values in edgeFilterProp
    g.set_edge_filter(edgeFilterProp)

    # Make network image using
    graph_draw(g, vertex_text_color='r', vertex_fill_color=drugSensProp3, vertex_size=10, vertex_text = cellLineProp, vertex_font_size=12, vertex_halo=True, vertex_halo_color = ccProp3, output_size=(4000,4000), output=figFolder + expTag + 'CellLineNetworkModel_STD.png')

    g.set_edge_filter(None)

    return

def GetSubspaceRelationship(g, highSensInds, lowSensInds):

    # Extract adjacency matrix and make sure that diagonals are = 1.0
    adj = adjacency(g, weight=g.edge_properties['simProp']).todense()
    for ii in range(g.num_vertices()):
        adj.itemset((ii,ii),1.0)

    # Collect highSens relvecs
    highSensVecs = list()
    # clinVals = g.vertex_properties['clinVals']
    for ii in highSensInds:
        curVec = list()
        for jj in range(len(adj)):
            curVec.append(adj.item((jj,ii)))
        highSensVecs.append(curVec)

    # Collect lowSens relvecs
    lowSensVecs = list()
    for ii in lowSensInds:
        if ii in lowSensInds:
            curVec = list()
            for jj in range(len(adj)):
                curVec.append(adj.item((jj,ii)))
            lowSensVecs.append(curVec)

    # Get a measure of the overlap of the two subspaces, weighted by the singular values
    uint,sint,vint = np.linalg.svd(np.dot(np.array(highSensVecs), np.array(lowSensVecs).transpose()))
    un, sn, vn = np.linalg.svd(np.dot(np.array(lowSensVecs),np.array(lowSensVecs).transpose()))

    # Sum up and take ratio of svals
    highSensFracPower = sum(sint) / sum(sn)

    return highSensFracPower

    # #
    # # DO DRUG-SPECIFIC ANALYSIS OF GRAPH NETWORK
    # #
    #
    # # Visualize graph
    # drug3SensProp = g.new_vertex_property('float')
    # for ii in range(g.num_vertices()):
    #     if drug3Sens[ii] != 'NA':
    #         drug3SensProp[g.vertex(ii)] = float(drug3Sens[ii])
    #     else:
    #         drug3SensProp[g.vertex(ii)] = 0.0
    #
    # sims = list()
    # for ii in range(g.num_vertices()):
    #     for jj in range(g.num_vertices()):
    #         sims.append(simProp[g.edge(g.vertex(ii),g.vertex(jj))])
    #
    # # Get binary clinical value as to whether mutation was present or not
    # clinValsRaw = g.vertex_properties['clinVals']
    # clinVals = g.new_vertex_property('float')
    # for ii in range(g.num_vertices()):
    #     # clinVals.append( clinValsRaw[g.vertex(ii)][0] )
    #     clinVals[g.vertex(ii)] = clinValsRaw[g.vertex(ii)][0]
    #
    # # Create edge filter to eliminate edges and turn off edges
    # edgeFilterProp = g.new_edge_property('bool')
    # for ii in range(g.num_vertices()):
    #     for jj in range(g.num_vertices()):
    #         if clinVals[g.vertex(ii)] == 1 and clinVals[g.vertex(jj)] == 1:
    #             edgeFilterProp[g.edge(g.vertex(ii),g.vertex(jj))]=simProp[g.edge(g.vertex(ii),g.vertex(jj))]
    # g.set_edge_filter(edgeFilterProp)
    #
    # # # Make edge histogram
    # # counts, bins = graph_tool.stats.edge_hist(g, simProp, [0,0.25,0.50,0.75,1.0,1.25,1.5])
    # # # plt.hist()
    #
    # # degTotal = g.degree_property_map('total')
    # # simProp = g.edge_properties['simProp']
    #
    # # Make vertex feature for drug 3 sens
    # drug3SensProp = g.new_vertex_property('float')
    # for ii in range(g.num_vertices()):
    #     if ii in highSens3:
    #         drug3SensProp[g.vertex(ii)] = 1
    #     elif ii in lowSens3:
    #         drug3SensProp[g.vertex(ii)] = 0
    #     else:
    #         drug3SensProp[g.vertex(ii)] = 0.5
    #
    # graph_draw(g, vertex_text_color='r', vertex_fill_color=clinVals, vertex_size=10, vertex_text = cellLineProp, vertex_font_size=12, vertex_halo=True, vertex_halo_color = drug3SensProp, output_size=(4000,4000), output=figFolder + 'CellLineNetworkModel_STD.png')
    #
    # # Create graph showing SVD-based segmentation of cell line using halo color
    # svdClassProp = g.new_vertex_property('bool')
    # for ii in range(g.num_vertices()):
    #     if u1[ii] >= 0:
    #         svdClassProp[g.vertex(ii)] = 1
    #     else:
    #         svdClassProp[g.vertex(ii)] = 0
    #
    # # g.set_edge_filter(None)
    # graph_draw(g, vertex_text_color='r', vertex_fill_color=clinVals, vertex_size=10, vertex_text = cellLineProp, vertex_font_size=12, vertex_halo=True, vertex_halo_color = svdClassProp, output_size=(4000,4000), output=figFolder + 'CellLineNetworkModel_STD_SVDClass.png')
    #
    # return


    #
    # # Get SVD of each
    # PrintStartTime('Started calculating SVDs ')
    # u0,s0,v0 = np.linalg.svd(adj0)
    # u1,s1,v1 = np.linalg.svd(adj1)
    #
    # # Use second singular vector to segment population
    # u20 = u0[:,1]
    # u21 = u0[:,1]
    #
    # ctr = 0
    # for ii in range(len(u20)):
    #     if u20[ii] > 0:
    #        ctr += 1
    #
    # # Get out clin val from vertices
    # clinValNeg = list()
    # clinValPos = list()
    # clinValProp = g0.vertex_properties['clinVals']
    # for ii in range(g0.num_vertices()):
    #     if u20[ii] > 0:
    #         clinValPos.append(clinValProp[g0.vertex(ii)].dataValueList)
    #     else:
    #         clinValNeg.append(clinValProp[g0.vertex(ii)].dataValueList)
    #
    # # Re order the rows of both matrices to try to align the patients in
    # # each group into a common notional order
    #
    # sum0 = adj0.sum(0, 'float')
    # sum1 = adj1.sum(0, 'float')
    #
    # sumAdj0 = list()
    # sumAdj1 = list()
    # for ii in range(numPatients):
    #     sumAdj0.append(sum0[0,ii])
    #     sumAdj1.append(sum1[0,ii])
    #
    # sortInds0 = sorted(range(len(sumAdj0)), key=lambda k: sumAdj0[k])
    # sortInds1 = sorted(range(len(sumAdj1)), key=lambda k: sumAdj1[k])
    #
    # e0 = np.zeros((numPatients, numPatients))
    # e1 = np.zeros((numPatients, numPatients))
    # for ii in range(numPatients):
    #     e0.itemset((ii,sortInds0[ii]), 1.0)
    #     e1.itemset((ii,sortInds1[ii]), 1.0)
    #
    # rotAdj0 = e0 * adj0 * e0.transpose()
    # rotAdj1 = e1 * adj1 * e1.transpose()
    #
    # uurot0, srot0, vrot0 = np.linalg.svd(rotAdj0)
    # uurot1, srot1, vrot1 = np.linalg.svd(rotAdj1)
    #
    # #

def ExtractGATally():
    """

        Method to determine whether there any unique or common GAs
        in any of the subpopulations

    :return:
    """

    print 'begin'

    # Load subpop results
    f = open(resultsFolder + expTag + 'SubpopResults', 'r')
    numVerts = cPickle.load(f)
    numVarsUsed = cPickle.load(f)
    gName = cPickle.load(f)
    varsUsed = cPickle.load(f)
    clinVals = cPickle.load(f)
    cellLines = cPickle.load(f)
    gListClosure = cPickle.load(f)
    f.close()

    gList, gListClosure = RebuildGList()
    gListInd = [ ii for ii,xx in enumerate(gListClosure) if xx == 0 ]

    # Load raw graph data to get list of total cell lines:  cellLinesX
    f1 = open(graphsFolder + expTag + 'RawGraphData','r')
    cellLinesX = cPickle.load(f1)
    drugSensValsX = cPickle.load(f1)
    valuesList = cPickle.load(f1)
    studyLineAbbX = cPickle.load(f1)
    f1.close()

    # Iterate thru each file and pull out values, mutations, and gene Name
    cellLineAltVal = list()
    cellLineAlt = list()
    geneName = list()
    for ii in range(len(geneAlterationFiles)):
        # Read in mutations information
        f0 = open(savedDataFolder + expTag + geneAlterationFiles[ii][1],'r')
        cellLineAltVal.append(cPickle.load(f0))
        cellLineAlt.append(cPickle.load(f0))
        geneName.append(sorted(cPickle.load(f0)))
        f0.close()

    #
    # PROFILE CELL LINES TO GET DATA ABOUT GENES THAT EACH LINES HAS MUT/CNVS FOR
    #
    cellLineActiveAlts, cellLineActiveVals, numAlts, geneAlts = ProfileCellLines(cellLinesX, cellLineAltVal, cellLineAlt, geneName)

    numSubpops = len(numVerts)

    # Iterate through each subpop and look at active GAs among cell lines in the
    # subpop
    subpopGas = list()
    subpopIds = list()
    subpopMutTally = list()
    subpopCnvTally = list()
    masterMutTally = Counter()
    masterCnvTally = Counter()
    for ii in range(numSubpops):

        curG = gList[gListInd[ii]]
        curIdProp = curG.vertex_properties['idProp']

        # Get cell line IDs
        curCellLines = list()
        curSubpopIds = list()
        for jj in range(curG.num_vertices()):
            curCellLines.append(curIdProp[curG.vertex(jj)])
            curSubpopIds.append(curIdProp[curG.vertex(jj)])

        # Find corresponding GAs
        curMuts = list()
        curCnvs = list()
        curCnvVals = list()
        curMutTally = Counter()
        curCnvTally = Counter()
        for jj in range(len(curCellLines)):
            # find index of current cell line
            curInd = [ kk for kk,xx in enumerate(cellLinesX) if curCellLines[jj] == xx ]
            curMuts.append(cellLineActiveAlts[curInd[0]][0])
            curCnvs.append(cellLineActiveAlts[curInd[0]][1])
            curCnvVals.append(cellLineActiveVals[curInd[0]][1])

            # Tally the occurrence of mutations and CNVs within each subpop
            for kk in range(len(curMuts[jj])):
                curMutTally[curMuts[jj][kk]] += 1
                masterMutTally[curMuts[jj][kk]] += 1
            for kk in range(len(curCnvs[jj])):
                curCnvTally[curCnvs[jj][kk]] += 1
                masterCnvTally[curCnvs[jj][kk]] += 1

        subpopMutTally.append(curMutTally)
        subpopCnvTally.append(curCnvTally)

        subpopGas.append( [ curMuts, curCnvs, curCnvVals ])
        subpopIds.append(curSubpopIds)

    # Iterate thru each MUT/CNV in masterCnvTally and compare counts to
    # fraction of counts in subpop Mut/Cnv tallycounts

    chiSqMaster = list()
    pValMaster = list()
    countsMaster = list()
    obsMaster = list()
    mutMaster = list(set(masterMutTally.elements()))
    for curMut in mutMaster:

        # Get total counts of current gene in masterMutTally
        curMasterCounts = float(masterMutTally[curMut])

        if curMasterCounts > 10:

            #
            expTemp = curMasterCounts * np.divide(np.array(numVerts), float(sum(numVerts)))
            curObs = list()
            for jj in range(len(numVerts)):
                curObs.append(float(subpopMutTally[jj][curMut]))

            chisq, p = chisquare(np.array(curObs), expTemp)
            chiSqMaster.append(chisq)
            pValMaster.append(p)
            countsMaster.append(curMasterCounts)
            obsMaster.append(curObs)

        else:
            chiSqMaster.append(0)
            pValMaster.append(1)
            countsMaster.append(curMasterCounts)
            obsMaster.append([])
    print 'end'
    f1 = open(resultsFolder + expTag + 'GeneticAlterationTallies','wb')
    cPickle.dump(subpopMutTally,f1)
    cPickle.dump(subpopCnvTally,f1)
    cPickle.dump(subpopGas,f1)
    cPickle.dump(subpopIds,f1)
    cPickle.dump(masterMutTally,f1)
    cPickle.dump(masterCnvTally,f1)
    f1.close()


    return

def ExtractGeneticAlterations():
    """

        Method to determine whether there any unique or common GAs
        in any of the subpopulations

    :return:
    """

    # Load subpop results
    f = open(localResultsFolder + expTag + 'SubpopResults', 'r')
    numVerts = cPickle.load(f)
    numVarsUsed = cPickle.load(f)
    gName = cPickle.load(f)
    varsUsed = cPickle.load(f)
    clinVals = cPickle.load(f)
    cellLines = cPickle.load(f)
    gListClosure = cPickle.load(f)
    f.close()

    gList, gListClosure = RebuildGList()
    gListInd = [ ii for ii,xx in enumerate(gListClosure) if xx == 0 ]

    # Load raw graph data to get list of total cell lines:  cellLinesX
    f1 = open(graphsFolder + expTag + 'RawGraphData','r')
    cellLinesX = cPickle.load(f1)
    drugSensValsX = cPickle.load(f1)
    valuesList = cPickle.load(f1)
    studyLineAbbX = cPickle.load(f1)
    f1.close()

    # Iterate thru each file and pull out values, mutations, and gene Name
    cellLineAltVal = list()
    cellLineAlt = list()
    geneName = list()
    for ii in range(len(geneAlterationFiles)):
        # Read in mutations information
        f0 = open(savedDataFolder + expTag + geneAlterationFiles[ii][1],'r')
        cellLineAltVal.append(cPickle.load(f0))
        cellLineAlt.append(cPickle.load(f0))
        geneName.append(sorted(cPickle.load(f0)))
        f0.close()

    #
    # PROFILE CELL LINES TO GET DATA ABOUT GENES THAT EACH LINES HAS MUT/CNVS FOR
    #
    cellLineActiveAlts, cellLineActiveVals, numAlts, geneAlts = ProfileCellLines(cellLinesX, cellLineAltVal, cellLineAlt, geneName)

    numSubpops = len(numVerts)

    # Iterate through each subpop and look at active GAs among cell lines in the
    # subpop
    subpopGas = list()
    subpopIds = list()
    for ii in range(numSubpops):

        curG = gList[gListInd[ii]]
        curIdProp = curG.vertex_properties['idProp']

        # Get cell line IDs
        curCellLines = list()
        curSubpopIds = list()
        for jj in range(curG.num_vertices()):
            curCellLines.append(curIdProp[curG.vertex(jj)])
            curSubpopIds.append(curIdProp[curG.vertex(jj)])

        # Find corresponding GAs
        curMuts = list()
        curCnvs = list()
        curCnvVals = list()
        for jj in range(len(curCellLines)):
            # find index of current cell line
            curInd = [ kk for kk,xx in enumerate(cellLinesX) if curCellLines[jj] == xx ]
            curMuts.append(cellLineActiveAlts[curInd[0]][0])
            curCnvs.append(cellLineActiveAlts[curInd[0]][1])
            curCnvVals.append(cellLineActiveVals[curInd[0]][1])
        subpopGas.append( [ curMuts, curCnvs, curCnvVals ])
        subpopIds.append(curSubpopIds)

    # Identify unique genetic drivers for each subgroup with respect to rest of
    # cell lines
    uniqGas = list()
    commGas = list()
    for ii in range(numSubpops):
        # Find GAs that are unique and distinct to a subpop compared to all other ones
        compMuts = list()
        compCnvs = list()
        compPopInds = range(numSubpops)
        compPopInds.pop(ii)
        for jj in compPopInds:
            compMuts = compMuts + subpopGas[jj][0]
            compCnvs = compCnvs + subpopGas[jj][1]

        # Make mutations from complementary set and current subpop into sets
        compMutSet = set(FlattenList(compMuts))
        subpopMutSet = set(FlattenList(subpopGas[ii][0]))

        # Make Cnvs from complementary set and current subpop into sets
        compCnvsSet = set(FlattenList(compCnvs))
        subpopCnvsSet = set(FlattenList(subpopGas[ii][1]))

        # Take difference to find out if anything in the subpop sets is unique to the complement
        # curUniqMuts = set(FlattenList(compMuts)).difference_update(set(FlattenList(subpopGas[ii][0])))
        # curUniqCnvs = set(FlattenList(compCnvs)).difference_update(set(FlattenList(subpopGas[ii][1])))
        uniqGas.append([subpopMutSet.difference_update(compMutSet), subpopCnvsSet.difference_update(compCnvsSet)])

        # Find the GAs that are common within each subpop
        # for jj in range(0,len(subpopGas[ii])):
        curCommMuts = subpopGas[ii][0][0]
        curCommCnvs = subpopGas[ii][1][0]
        curCommCnvVals = subpopGas[ii][2][0]
        x=5
        # Iterate thru each subsequent cell line in the current subgroup
        for kk in range(1,len(subpopGas[ii][0])):
            # print kk
            # Get xsection of mutations
            curCommMuts = list(set(curCommMuts).intersection(set(subpopGas[ii][0][kk])))

            # Get common CNVS based on name
            # masterCommCnvList = list(curCommCnvs)
            for nn in reversed(range(len(curCommCnvs))):
                if curCommCnvs[nn] not in subpopGas[ii][1][kk]:
                    curCommCnvs.pop(nn)
                    curCommCnvVals.pop(nn)

            # curCommCnvs = list(set(curCommCnvs).intersection(set(subpopGas[ii][1][kk])))
            # for mm in reversed(range(len(masterCommCnvList))):
            #     if masterCommCnvList[mm] not in curCommCnvs:
            #         curCommCnvVals.pop(mm)

            # Filter remaining CNVs by value
            # print '**', ii, kk, len(curCommCnvVals), len(curCommCnvs)
            popList = list()
            for mm in range(len(curCommCnvs)):
                # Get index of cnv value for common set
                curCommCnvInd = [ zz for zz,xx in enumerate(curCommCnvs) if curCommCnvs[mm] == xx  ][0]
                curSubpopCnvInd = [ zz for zz,xx in enumerate(subpopGas[ii][1][kk]) if curCommCnvs[mm] == xx  ][0]
                # print mm, curCommCnvInd, curSubpopCnvInd, len(curCommCnvVals), len(subpopGas[ii][2][kk])
                if curCommCnvVals[curCommCnvInd] != subpopGas[ii][2][kk][curSubpopCnvInd]:
                    popList.append(curCommCnvInd)
                    # print mm, curCommCnvs[mm], curCommCnvVals[curCommCnvInd], subpopGas[ii][2][kk][curSubpopCnvInd]

            for mm in sorted(popList, reverse=True):
                curCommCnvVals.pop(mm)
                curCommCnvs.pop(mm)

        commGas.append([curCommMuts, curCommCnvs, curCommCnvVals ])

    f0 = open(resultsFolder + expTag + 'DistinctGeneticAlterations','wb')
    cPickle.dump(uniqGas, f0)
    cPickle.dump(commGas, f0)
    f0.close()

    return uniqGas, commGas

def MakeSubpopBoxplots():
    """

        Method to make a boxplot summary of all the subpops that were ID'd
        in terms of their AUC and IC50 performance

    :return:
    """

    # Retrieve data
    # gList, gListClosure = RebuildGList()
    # numVerts, numVarsUsed, gName, varsUsed, clinVals, cellLines, gListInd = ExtractSubpopulations(gList, gListClosure)

    # Write results to file
    f = open(resultsFolder + expTag + curTarg + '_' + 'SubpopResults', 'r')
    numVerts = cPickle.load(f)
    numVarsUsed = cPickle.load(f)
    gName = cPickle.load(f)
    varsUsed = cPickle.load(f)
    clinVals = cPickle.load(f)
    cellLines = cPickle.load(f)
    gListClosure = cPickle.load(f)
    f.close()

    numSubpops = len(numVerts)
    numCovs = len(varsUsed[0])

    if numSubpops <= 8:

        # boxplot of clustering results
        fig = plt.figure(figsize=(14,10))

        # Determine number of rows in subplot page
        numRows = math.ceil(np.divide(float(numSubpops),4.0))

        # Get min and max value for boxplots
        # Turn clinVals into one long list
        clinValsList = list()
        for ii in range(len(clinVals)):
            for jj in range(len(clinVals[ii])):
                for kk in range(len(clinVals[ii][jj])):
                    clinValsList.append(clinVals[ii][jj][kk])

        minY = min(clinValsList)
        maxY = max(clinValsList)

        # Make list of string labels
        respVarLabels = list()
        for ii in range(len(allVarsInclLabels)):
            respVarLabels.append(allVarsInclLabels[ii].replace('_',' \n '))

        # Get membership of each tissue type in each subpop
        if tissType[0] == 'ALL':
            tissTypeMem = list()
            for ii in range(numSubpops):
                tissTypeMem.append('ALL')
        else:
            tissTypeMem = list()
            for ii in range(numSubpops):
                tissTypeCtr = [0] * len(tissType)
                curG = gList[gListInd[ii]]
                studyLineAbbProp = curG.vertex_properties['studyLineAbb']
                for jj in range(gList[gListInd[ii]].num_vertices()):
                    curStudyLine = studyLineAbbProp[curG.vertex(jj)]
                    tissTypeInd = [ kk for kk,xx in enumerate(tissType) if xx == curStudyLine ]
                    # tissTypeCtr[tissTypeInd[0]] += 1
                    tissTypeCtr[tissTypeInd[0]] += 1
                tissTypeMem.append(tissTypeCtr)

        # Produce boxplot using resulting clinVals from each subpop
        ax = list()
        for ii in range(len(clinVals)):

            ax.append(fig.add_subplot(numRows,4,ii+1))
            plt.boxplot(PivotList(clinVals[ii]))
            ax[ii].set_xticks(list(np.arange(1.0,numCovs+2,1)))
            ax[ii].set_xticklabels(allVarsInclLabels, rotation='vertical', horizontalalignment='center')

            # ax[ii].set_xlabel('Phenotypic Drug Response Variables')
            ax[ii].set_ylabel('Output Values (AUC, log(IC50)')
            x1,x2,y1,y2 = ax[ii].axis()
            ax[ii].axis([x1,x2,minY,maxY])

            plt.grid()
            ax[ii].set_title('Subpop ' + str(ii+1) + ' (N = ' + str(numVerts[ii]) + ')\n' + 'Tissue Memb. = ' + str(tissTypeMem[ii]))
        plt.tight_layout()

        # Make list of string labels
        tissTypeLabel = tissType[0]
        for ii in range(1,len(tissType)):
            tissTypeLabel = tissTypeLabel + '_' + tissType[ii]

        plt.savefig(figFolder + expTag + 'CellLineResponseClustering_' + tissTypeLabel)
        plt.close()

        # Plot mean AUC and IC50 values for AZD 9291 for each subpop
        meanAuc = list()
        meanIc50 = list()
        stdAuc = list()
        stdIc50 = list()
        for ii in range(len(clinVals)):
            meanAuc.append(np.mean(PivotList(clinVals[ii])[0]))
            meanIc50.append(np.mean(PivotList(clinVals[ii])[1]))
            stdAuc.append(np.std(PivotList(clinVals[ii])[0]))
            stdIc50.append(np.std(PivotList(clinVals[ii])[1]))

        fig = plt.figure()
        plt.bar(np.arange(0.6,numSubpops+0.5,1), meanAuc, 0.4, color='b', yerr=stdAuc, label=drugName[0] + ' Mean AUC')
        plt.bar(np.arange(1,numSubpops+1,1), meanIc50, 0.4, bottom=0, color='r', yerr=stdIc50, label=drugName[0] + ' Mean IC50')
        plt.grid()
        plt.xlabel('Subpop. No.')
        plt.ylabel('Mean Response Value')
        plt.title('Mean Drug Response for ' + drugName[0] + ' for Subpops. \n Studies included are ' + str(tissType))
        plt.legend(loc=0)
        plt.tight_layout()
        plt.savefig(localResultsFolder + expTag + curTarg + '_' + drugName[0] + ' Mean Subpop Response.png')
        plt.close()

    else:
        # Get min and max value for boxplots
        # Turn clinVals into one long list
        clinValsList = list()
        for ii in range(len(clinVals)):
            for jj in range(len(clinVals[ii])):
                for kk in range(len(clinVals[ii][jj])):
                    clinValsList.append(clinVals[ii][jj][kk])

        minY = min(clinValsList)
        maxY = max(clinValsList)

        # Make list of string labels
        respVarLabels = list()
        for ii in range(len(allVarsInclLabels)):
            respVarLabels.append(allVarsInclLabels[ii].replace('_',' \n '))

        # Get membership of each tissue type in each subpop
        if tissType[0] == 'ALL':
            tissTypeMem = list()
            for ii in range(numSubpops):
                tissTypeMem.append('ALL')
        else:
            tissTypeMem = list()
            for ii in range(numSubpops):
                tissTypeCtr = [0] * len(tissType)
                curG = gList[gListInd[ii]]
                studyLineAbbProp = curG.vertex_properties['studyLineAbb']
                for jj in range(gList[gListInd[ii]].num_vertices()):
                    curStudyLine = studyLineAbbProp[curG.vertex(jj)]
                    tissTypeInd = [ kk for kk,xx in enumerate(tissType) if xx == curStudyLine ]
                    # tissTypeCtr[tissTypeInd[0]] += 1
                    tissTypeCtr[tissTypeInd[0]] += 1
                tissTypeMem.append(tissTypeCtr)

        # Produce boxplot using resulting clinVals from each subpop
        for ii in range(numSubpops):

            # boxplot of clustering results
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)

            plt.boxplot(PivotList(clinVals[ii]))
            ax.set_xticks(list(np.arange(1.0,numCovs+2,1)))
            ax.set_xticklabels(allVarsInclLabels, rotation='vertical', horizontalalignment='center')

            # ax[ii].set_xlabel('Phenotypic Drug Response Variables')
            ax.set_ylabel('Output Values (AUC, log(IC50)')
            x1,x2,y1,y2 = ax.axis()
            ax.axis([x1,x2,minY,maxY])

            plt.grid()
            ax.set_title('Subpop ' + str(ii+1) + ' (N = ' + str(numVerts[ii]) + ')\n' + 'Tissue Memb. = ' + str(tissTypeMem[ii]))
            plt.tight_layout()

            # Make list of string labels
            tissTypeLabel = tissType[0]
            for ii in range(1,len(tissType)):
                tissTypeLabel = tissTypeLabel + '_' + tissType[ii]

            plt.savefig(localResultsFolder + expTag + curTarg + '_' + 'CellLineResponseCluster_' + str(ii+1) + '_' + tissTypeLabel)
            plt.close()

        # Plot mean AUC and IC50 values for AZD 9291 for each subpop
        meanAuc = list()
        meanIc50 = list()
        stdAuc = list()
        stdIc50 = list()
        for ii in range(len(clinVals)):
            meanAuc.append(np.mean(PivotList(clinVals[ii])[0]))
            meanIc50.append(np.mean(PivotList(clinVals[ii])[1]))
            stdAuc.append(np.std(PivotList(clinVals[ii])[0]))
            stdIc50.append(np.std(PivotList(clinVals[ii])[1]))

        fig = plt.figure()
        plt.bar(np.arange(0.6,numSubpops+0.5,1), meanAuc, 0.4, color='b', yerr=stdAuc, label=drugName[0] + ' Mean AUC')
        plt.bar(np.arange(1,numSubpops+1,1), meanIc50, 0.4, bottom=0, color='r', yerr=stdIc50, label=drugName[0] + ' Mean IC50')
        plt.grid()
        plt.xlabel('Subpop. No.')
        plt.ylabel('Mean Response Value')
        plt.title('Mean Drug Response for ' + drugName[0] + ' for Subpops. \n Studies included are ' + str(tissType))
        plt.legend(loc=0)
        plt.tight_layout()
        plt.savefig(localResultsFolder + expTag + curTarg + '_' + drugName[0] + ' Mean Subpop Response.png')
        plt.close()

        # Extract  logIc50 data for the 3 drugs, in each subpop
        meanIc50 = [ [] for ii,xx in enumerate(range(numDrugs)) ]
        stdIc50 = [ [] for ii,xx in enumerate(range(numDrugs)) ]
        for ii in range(len(clinVals)):
            for jj in range(numDrugs):
                meanIc50[jj].append(np.mean(PivotList(clinVals[ii])[2*jj + 1]))
                # meanIc50[1].append(np.mean(PivotList(clinVals[ii])[3]))
                # meanIc50[2].append(np.mean(PivotList(clinVals[ii])[5]))
                stdIc50[jj].append(np.std(PivotList(clinVals[ii])[2*jj + 1]))
                # stdIc50[1].append(np.std(PivotList(clinVals[ii])[3]))
                # stdIc50[2].append(np.std(PivotList(clinVals[ii])[5]))

        MakeIc50Subpop(meanIc50, clinVals, numVerts)

    return

def QuantizeSubpops():

    # Write results to file
    f = open(localResultsFolder + expTag + curTarg + '_' + 'SubpopResults', 'r')
    numVerts = cPickle.load(f)
    numVarsUsed = cPickle.load(f)
    gName = cPickle.load(f)
    varsUsed = cPickle.load(f)
    clinVals = cPickle.load(f)
    cellLines = cPickle.load(f)
    gListClosure = cPickle.load(f)
    f.close()
    Nsubpops = len(numVerts)

    # Apply

    # Load data from composite variable analysis; most importantly
    # retrieve the indexed order of most sensitive --> most resistant drug
    f = open(resultsFolder + expTag + curTarg + '_CompositeSensitivityVariables', 'r')
    compVals = cPickle.load(f)
    sortSubpopInds = cPickle.load(f)
    f.close()

    # Get mean comp val for each drug
    meanCompVal = list()
    medCompVal = list()
    stdCompVal = list()
    for ii in range(len(compVals)):
        meanCompVal.append([])
        medCompVal.append([])
        stdCompVal.append([])
        for jj in range(len(compVals[ii])):
            meanCompVal[ii].append(np.mean(compVals[ii][jj]))
            medCompVal[ii].append(np.median(compVals[ii][jj]))
            stdCompVal[ii].append(np.std(compVals[ii][jj]))

    plt.scatter(meanCompVal[0], meanCompVal[1], color='b', label='Mean')
    plt.scatter(medCompVal[0], medCompVal[1], color='r', label='Median')
    plt.xlabel(drugName[0])
    plt.ylabel(drugName[1])
    plt.title('Scatterplot of Composite Variable for ' + str(Nsubpops) + ' Subpops for \n ' + str(drugName))
    plt.grid()
    plt.legend(loc='best')
    plt.savefig(figFolder + expTag + 'MeanMedianCompValScatterPlot.png')
    plt.close()


    # Q subpops into four groups


    x = 5


    return

def ProfileCellLines(cellLineKeep, cellLineAltVal, cellLineAlt, geneName):
    """

        Method to profile a list of cell lines and return information on which genetic mutations
        are present in each cell line and how many are presesnt in the aggregate

    :param cellLineKeep:
    :param geneAlterationFiles:
    :return:
    """

    # Iterate thru each cell line in cellLineKeep and collect all genes and mutations
    # for each cell line
    cellLineActiveAlts = list()
    cellLineActiveVals = list()
    for ii in range(len(cellLineKeep)):
        # Get current cell line
        curCellLine = cellLineKeep[ii]

        # Find index of entry in cellLineAlt entries and collect all gene values for
        # each cell line
        curCellLineInd = list()
        curGeneVals = list()
        for jj in range(len(cellLineAlt)):
            curInd = [kk for kk,xx in enumerate(cellLineAlt[jj]) if xx == curCellLine]
            curCellLineInd.append(curInd[0])
            curGeneVals.append(cellLineAltVal[jj][curCellLineInd[jj]])

        # Get total list of non-zero mutations in each cell line of all possible types
        curActiveAltsList = list()
        curActiveValsList = list()
        for jj in range(len(cellLineAlt)):
            # Get the list of genes in the current alt file that have non zero values
            curActiveAltsList.append([ geneName[jj][kk] for kk,xx in enumerate(cellLineAltVal[jj][curCellLineInd[jj]]) if xx != 0 ])
            # Get the list of the value of genes in the current alt file that have non zero values
            curActiveValsList.append([ xx for kk,xx in enumerate(cellLineAltVal[jj][curCellLineInd[jj]]) if xx != 0 ])

        cellLineActiveAlts.append(curActiveAltsList)
        cellLineActiveVals.append(curActiveValsList)

    #
    # GENERATE STATISTICS ABOUT GENETIC ALTERATIONS IN CELL LINES
    #

    # Number of alterations of each type, by cell line
    numAlts = list()
    for ii in range(len(cellLineKeep)):
        curAlts = cellLineActiveAlts[ii]
        curLens = list()
        for jj in range(len(curAlts)):
            curLens.append(len(curAlts[jj]))
        numAlts.append(curLens)

    # Distinct list of genes active in each alteration type, across entire cell line
    geneAlts = list()
    for ii in range(len(numAlts[0])):
        geneAlts.append([])
    for ii in range(len(cellLineActiveAlts)):
        for jj in range(len(cellLineActiveAlts[ii])):
            for kk in range(len(cellLineActiveAlts[ii][jj])):
                geneAlts[jj].append(cellLineActiveAlts[ii][jj][kk])

    for ii in range(len(geneAlts)):
        geneAlts[ii] = list(set(geneAlts[ii]))

    return cellLineActiveAlts, cellLineActiveVals, numAlts, geneAlts

def PartitionGraphAndSegment(g, parentLabel):
    """

        Method to partition a graph, g, having label, parentLabel
        into 2 graphs using the Fiedler eigenvector

    :param g:
    :param parentLabel:
    :return:
    """

    # print'***************************************************************************'
    # PrintStartTime('Currently partitioning graph, ' + parentLabel + ' having ' + str(g.num_vertices()) + ' vertices.')

    #
    # DERIVE FIEDLER EIGENVECTOR AND SEGMENT CLN
    #
    adj = adjacency(g, g.edge_properties['simProp']).todense()
    lap = PerformLaplacian(adj)

    # Perform svd of lap
    u,s,v = np.linalg.svd(lap)
    fVec = u[:,len(s)-2]
    fVal = s[len(s)-2]

    # Convert fVec into binary vector where vals < 0 --> 0 and vals > 0 --> 1
    binFVec = list()
    for jj in range(len(fVec)):
        if fVec[jj] > 0.0:
            binFVec.append(1)
        else:
            binFVec.append(0)
    if sum(binFVec) > math.floor(len(binFVec)/2.0):
        for ii in range(len(binFVec)):
            if binFVec[ii] == 0:
                binFVec[ii] = 1
            else:
                binFVec[ii] = 0

    # Extract clinical values for class0 and class1
    clinValsProp = g.vertex_properties['clinVals']
    clinVals = ExtractListsFromVertices(clinValsProp, g)
    # Segregate clinVals into two groups, class0 and class1
    clinVals0 = [ clinVals[ii] for ii,xx in enumerate(binFVec) if xx == 0 ]
    clinVals1 = [ clinVals[ii] for ii,xx in enumerate(binFVec) if xx == 1 ]

    # Extract id property and segregate
    idProp = g.vertex_properties['idProp']
    ids = ExtractListsFromVertices(idProp,g)
    ids0 = [ ids[ii] for ii,xx in enumerate(binFVec) if xx == 0 ]
    ids1 = [ ids[ii] for ii,xx in enumerate(binFVec) if xx == 1 ]

    # print 'Fiedler eigenvalue is ' + str(fVal)
    # print 'Partitioning found ' + str(len(ids0)) + ' vertices in first partition of ' + parentLabel
    # print 'Partitioning found ' + str(len(ids1)) + ' vertices in second partition of ' + parentLabel

    #
    # Calculate silhouette score for partition
    #

    # Get adjacency from g
    adj = np.array(adjacency(g, g.edge_properties['simProp']).todense())
    # Make sure diagonal values are = 1.0
    for ii in range(g.num_vertices()):
        adj.itemset((ii,ii), 1.0)
    sscore = silScore(1.0 - adj, np.array(binFVec), metric='precomputed')
    # print 'Silhouette score of new graphs is ' + str(sscore)

    # Create new graphs for class0 and class1
    g0, curMatsList0 = MakeMultivariatePatientNetworkGraph(PivotList(clinVals0), ids0)
    g1, curMatsList1 = MakeMultivariatePatientNetworkGraph(PivotList(clinVals1), ids1)

    # Get number ov variables actually used in graph construction
    varsUsed0 = g0.graph_properties['varsUsed']
    numVarsUsed0 = sum(varsUsed0)
    varsUsed1 = g1.graph_properties['varsUsed']
    numVarsUsed1 = sum(varsUsed1)

    # Add graph name to each graph
    graphName0 = g0.new_graph_property('string', parentLabel + '0')
    graphName1 = g1.new_graph_property('string', parentLabel + '1')

    # Add parent name to each graph
    parentName0 = g0.new_graph_property('string', parentLabel)
    parentName1 = g1.new_graph_property('string', parentLabel)

    print
    print 'Graph ' + graphName0[g0] + ' has ' + str(len(ids0)) + ' vertices using ' + str(numVarsUsed0) + ' variables'
    print 'Graph ' + graphName1[g1] + ' has ' + str(len(ids1)) + ' vertices using ' + str(numVarsUsed1) + ' variables'
    print
    print'***************************************************************************'
    print

    return g0, g1, binFVec, sscore, varsUsed0, varsUsed1

def GetDrugSensVals(cellLineIn, cellLinePh, drugSensVals, th):
    """

        Method to return the cell lines for which there are sensible values

    :param cellLineIn:
    :param cellLinePh:
    :param drugSensVals:
    :param th:
    :return:
    """

    sensVals = list()
    binSensVals = list()
    cellLineKeep = list()
    for ii in range(len(cellLineIn)):
        if cellLineIn[ii] in cellLinePh:
            # Find out index for current value of CellLineIn
            curInd = [jj for jj,xx in enumerate(cellLinePh) if cellLineIn[ii]==xx]
            if drugSensVals[curInd[0]] != 'NA':
                sensVals.append(float(drugSensVals[curInd[0]]))
                cellLineKeep.append(cellLineIn[ii])
                if float(drugSensVals[curInd[0]]) <= th:
                    binSensVals.append(1)
                else:
                    binSensVals.append(0)

    return cellLineKeep, sensVals, binSensVals

def ExtractVarFromTupleList(tupleList, n):
    """
        Method to extract from a list of tuples, one variable in the tuple and return the
        values in a list

    :param tupleList:
    :param n:
    :return:
    """

    varList = list()
    N = len(tupleList)

    for ii in range(N):
        varList.append(tupleList[ii][n])

    return varList

def ConvertListsToArray(varList):
    """
        Method to convert a list of lists (all of the same length) to one 2-D array

    :param varList:
    :return:
    """

    # Get dimensions of array from varList
    numVars = len(varList)
    N = len(varList[0])

    # Iterate thru each list and add populate an array
    varArray = np.ndarray((N, numVars))
    for jj in range(numVars):
        for kk in range(N):
            varArray[kk,jj] = varList[jj][kk]

    return varArray

def FlattenList(inputList):
    """

        Method to take a list of lists and produce one overall list

    :return:
    """

    opList = list()
    for ii in range(len(inputList)):
        opList = opList + inputList[ii]

    return opList
