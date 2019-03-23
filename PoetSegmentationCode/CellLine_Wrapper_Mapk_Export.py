__author__ = 'ktch802'

from CellLine_Wrapper_Mapk_Infofile import *

import time
import cPickle
import os
import math
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chisquare
from scipy.stats import ttest_ind as ttest
import csv
from graph_tool.all import *
from sklearn.metrics import silhouette_score as silScore

########################## UTILITIES #################

startTime = time.time()

def PrintStartTime(printString):
    """
        Print the current time since execution and a msg (printString)
    """

    endTime = time.time() - startTime
    print str(int(endTime)) + ' s: ' + printString + ' .......'

    return endTime

def PrintEndTime(printString):
    """
        Print the current time since execution and a msg (printString)
    """
    endTime = time.time() - startTime
    print str(int(endTime)) + ' s: ' + printString
    print

    return endTime

def MakeFolder(newFolder):
    """
        Method to make a new folder given the path of the desired folder;
        if the folder already exists, then it does nothing

    :param newFolder:
    :return:
    """

    if newFolder[-1] != '/':
        newFolder = newFolder + '/'

    # Make output folders if they don't exist
    if os.path.exists(newFolder):
        pass
    elif ~os.path.exists(newFolder):
        os.mkdir(newFolder)
        os.mkdir(newFolder + 'Figs/')
        os.mkdir(newFolder + 'SavedData/')
        os.mkdir(newFolder + 'Results/')
        os.mkdir(newFolder + 'Exports/')
        os.mkdir(graphsFolder)
        print 'Made ' + newFolder + ' and sub-folders'

    return newFolder

def PivotList(x):
    """
        Method to exchange dimensions of list of lists:  eg, list containing different
        variable values can be turned into list of patient lists, where each list contains
        multiple variables

    :param x:
    :return:
    """

    dim0 = len(x)

    if type(x[0]) == type([1]):
        newList = list()
        dim1 = len(x[0])
        for ii in range(dim1):
            curList = list()
            for jj in range(dim0):
                curList.append(x[jj][ii])
            newList.append(curList)

    else:
        newList = list()
        for ii in range(dim0):
            newList.append([x[ii]])

    return newList

def DeListAList(inputList):
    """
        Method to take a list in which each element is a list of length 1 and turn it into
        a list of elements

    :param inputList:  input list whose elements are lists of length-1
    :return:           opList
    """

    opList = list()
    for ii in range(len(inputList)):
        opList.append(inputList[ii][0])

    return opList

########################## READ IN CELL LINE DATA FROM XLS #################

def ReadInCellLineData(filename, opFilename):
    """

        Method to read in a .csv file and write it out to a python file and return its values

    :param filename:
    :param opFilename:
    :return:
    """

    # Open up csv file

    mutVal = list()
    gene = list()
    cellLine = list()
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        rowCtr = 0
        for row in reader:
            if rowCtr == 0:
                for kk in range(1,len(row)):
                    cellLine.append(row[kk])
            else:
                # Process each row to extract gene vs. value
                curRowVals = list()
                for jj in range(len(row)):
                    if jj == 0:
                        gene.append(row[0])
                    else:
                        curRowVals.append(int(row[jj]))
                mutVal.append(curRowVals)
                # print
            rowCtr += 1

    # Save results to a file, but transpose mutVals to be organized by cell line
    # (i.e., by patient)

    # Read in joint covariates and patient IDs
    f0 = open(savedDataFolder + opFilename,'wb')
    cPickle.dump(PivotList(mutVal), f0)
    cPickle.dump(cellLine, f0)
    cPickle.dump(gene, f0)
    f0.close()

    return PivotList(mutVal), gene, cellLine

def ReadPharmacologyData_AZD9291(filename, opFilename, numDrugs):

    # Open up csv file
    cellLine = list()
    # tissue = list()
    header = list()
    drugSens = list()
    for ii in range(numDrugs):
        drugSens.append([])
    # drug2Sens = list()
    # drug3Sens = list()
    sensData = list()
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        rowCtr = 0
        for row in reader:
            if rowCtr == 0:
                # Digest header row
                for kk in range(0,len(row)):
                    header.append(row[kk])

                # Find the column corresonding to entries in drugName
                drugNameColInd = list()
                for jj in range(len(drugName)):
                    drugNameColInd.append(header.index(drugName[jj]))

                rowCtr += 1
            else:
                # Capture rows as list
                curRow = list()
                for kk in range(0, len(row)):
                    curRow.append(row[kk])

                sensData.append(curRow)

                # Capture cell line in its own list
                # 5-5-16:  Do conversion of cell lines names to
                # remove hyphens
                cellLine.append(row[0].replace('-','').upper())
                # cellLine.append(row[0])
                # tissue.append(row[1])

                rowCtr += 1

    # Extract values for drugs in drugNames
    drugSens = list()
    for ii in range(len(drugName)):
        curVals = list()
        for jj in range(len(sensData)):
            curVals.append(sensData[jj][drugNameColInd[ii]])
        drugSens.append(curVals)

    # Save each drug response to a separate file with the same metadata
    drugSensFileNames = list()
    for ii in range(numDrugs):
        drugSensFileNames.append(savedDataFolder + expTag + opFilename + '_' + drugName[ii])
        f0 = open(drugSensFileNames[ii],'wb')
        cPickle.dump(header[drugNameColInd[ii]],f0)
        cPickle.dump(cellLine, f0)
        # cPickle.dump(tissue, f0)
        cPickle.dump(drugSens[ii], f0)
        f0.close()

    # return header, cellLine, tissue, drugSens, drugSensFileNames
    return header, cellLine, drugSens, drugSensFileNames

def ReadCellLineInfoFile_AZD9291(filename):

    # Open up csv file
    cellLine = list()
    cosmicId = list()
    studyAbb = list()
    cancerType = list()
    drugSens = list()
    header = list()
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        rowCtr = 0
        for row in reader:
            if rowCtr == 0:
                # Digest header row
                for kk in range(0,len(row)):
                    header.append(row[kk])
                rowCtr += 1
            else:
                # Capture rows as list
                curRow = list()
                for kk in range(0, len(row)):
                    curRow.append(row[kk])

                cellLine.append(curRow[0].replace('-','').upper())
                cosmicId.append(curRow[1])
                studyAbb.append(curRow[2])
                # cancerType.append(curRow[3])

                rowCtr += 1

    return cellLine, cosmicId, studyAbb, cancerType

def ReadCellLineInfoFile_Mapk(filename):

    # Open up csv file
    cellLine = list()
    cosmicId = list()
    studyAbb = list()
    cancerType = list()
    drugSens = list()
    header = list()
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        rowCtr = 0
        for row in reader:
            if rowCtr == 0:
                # Digest header row
                for kk in range(0,len(row)):
                    header.append(row[kk])
                rowCtr += 1
            else:
                # Capture rows as list
                curRow = list()
                for kk in range(0, len(row)):
                    curRow.append(row[kk])

                cellLine.append(curRow[1].replace('-','').upper())
                cosmicId.append(curRow[2])
                studyAbb.append(curRow[3])
                cancerType.append(curRow[4])

                rowCtr += 1

    return cellLine, cosmicId, studyAbb, cancerType

################### WRAPPER METHODS TO ENABLE LOOPING OVER MULTIPLE DRUG PAIRS #################

def ReadTargetsFile():
    """

        Method to read in the target list file that contiains for each
        entry the target of interest and the 2 compounds being considered for that
        target

    :return:
    """

    # Populate 2 lists:  targList holds list of targets
    #                    compList holds list of compounds associated with drugs
    #                    in targList
    targList = list()
    compList = list()
    with open(wrapperFile[0], 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            targList.append(row[0])
            curCompList = list()
            for kk in range(1,len(row)):
                if row[kk] != '':
                    curCompList.append(row[kk])
            compList.append(curCompList)

    return targList, compList

def IterateOverTargets(targList, compList, Ntargets):
    """

        Method to iterate over each tarets in the target list
        and execute POET for the associated compounds and cell lines

    :param targList:
    :param compList:
    :param Ntargets:
    :return:
    """

    # SEt up main loop with global variables that will change with every entry in the
    # wrapperFile
    global curTarg              # current Target under consideration
    global numDrugs             # number of drugs that are target agents for curTarg
    global drugName             # list of drug names that are agents for curTarg
    global allVarsInclLabels    # list of labels for variables being clustered
    global localResultsFolder   # local folder, based on curTarg, where output is
                                # directed

    # for ii in range(40,114):
    # for ii in range(114, 228):
    # for ii in range(228, 342):

    for ii in range(len(targList)):
        # print ii

        # Use ii-th target
        curTarg = targList[ii]

        # Assign global variables for current run of POET
        numDrugs = len(compList[ii])
        allVarsInclLabels = list()
        drugName = list()

        # Iterate thru each compound associatd with current target
        # to get drug name and assign string labels for each variable
        # Note:  drugSensitivityType refers to the list of variables
        #        that exist to describe the measured response of a cell line
        #        to a compound; examples are AUC and log(IC50)
        for jj in range(numDrugs):
            drugName.append(compList[ii][jj])
            for kk in range(numGeneAlterationTypes):
                allVarsInclLabels.append(drugName[jj] + drugSensitivityType[kk])

        # Define folder for current target and see if it exists; if not make it
        localResultsFolder = resultsFolder + curTarg + '/'
        if os.path.exists(localResultsFolder):
            pass
        elif ~os.path.exists(localResultsFolder):
            os.mkdir(localResultsFolder)

        ###########################################################################
        # Execute CellLine Master Method for current target and associated compounds
        # CellLine_Master_Method() is the method that creates a graph and partitions
        # it recursively
        ###########################################################################
        PrintStartTime(str(ii) + ' ; \t Starting POET for ' + curTarg + ' with ' + str(numDrugs) + ' compounds:' + str(compList[ii]))
        curNumVerts = CellLine_Master_Method()
        curNumVerts=[0]

        ###########################################################################
        # Print results of POET for current target
        ###########################################################################
        PrintEndTime(str(ii) + ' ; \t' + curTarg + ' ; \t' + str(len(drugName)) + ' drugs; \t' + str(sum(curNumVerts)) + ' cell lines \t' + str(len(curNumVerts)) + ' subpops ; \t')

        # DrawHierarchicalTree(figFolder + expTag + 'HierarchicalTree_' + str(ii+1) + '.png', resultsFolder+str(ii+1)+'/'+expTag+str(ii+1)+'_')
        DrawHierarchicalTreeWithOutcomes(figFolder + expTag + 'HierarchicalTree_' + str(ii+1) + '.png', resultsFolder+str(ii+1)+'/'+expTag+str(ii+1)+'_', '')

    return

############################## NETWORK GRAPH CONSTRUCTION ###############################

def MakeCellLineNetwork():
    """

        Method to create a cell line network (CLN) using the raw graph data values that have been
        processed, curated, and stored beforehand

    :param drugNum:
    :return:
    """

    # Save final lists of data for current processing
    f1 = open(graphsFolder + expTag + curTarg + '_' + 'RawGraphData','r')
    # f1 = open(graphsFolder + expTag + 'RawGraphData','r')
    cellLineX = cPickle.load(f1)
    sensValsX = cPickle.load(f1)
    # binSensVals = cPickle.load(f1)
    valuesList = cPickle.load(f1)
    studyLineAbbX = cPickle.load(f1)
    f1.close()

    # Make one list of lists with variables and then pivot
    sensValsList = list()
    for ii in range(len(drugName)):
        for jj in range(len(drugSensitivityType)):
            curVals = list()
            for kk in range(len(cellLineX)):
                curVals.append(sensValsX[ii][kk][jj])
            sensValsList.append(curVals)

    # Set graph-making variables to all floats
    global allVarsInclType
    allVarsInclType = list()
    global allVarsInclCats
    allVarsInclCats = list()
    for ii in range(len(drugName)*len(drugSensitivityType)):
        allVarsInclType.append('continuous')
        allVarsInclCats.append([])
    global allVarsInclLabels
    allVarsInclLabels = list()
    for ii in range(len(drugName)):
        for jj in range(len(drugSensitivityType)):
            allVarsInclLabels.append(drugName[ii]+'_'+drugSensitivityType[jj])

    #
    # CONSTRUCT GRAPH
    #

    # K=7
    # PerformKMeans(sensValsList, cellLineX, K)

    g0, curMatsList = MakeMultivariatePatientNetworkGraph(sensValsList, cellLineX)

    # Add studyLineAbbX as vertex property
    studyLineAbbXProp = g0.new_vertex_property('string')
    for ii in range(g0.num_vertices()):
        studyLineAbbXProp[g0.vertex(ii)] = studyLineAbbX[ii]
    g0.vertex_properties['studyLineAbb'] = studyLineAbbXProp

    # Save graph to file
    g0.save(graphsFolder + expTag + curTarg + '_' + 'CellLineNetwork.gt')

    return g0

def MakeMultivariatePatientNetworkGraph(valuesList, patientIdList):
    """
        Method to create multi-patient network model with measures of similarity
        between patients established using frequency counts of unique conditions
        that are shared pairwise between patients

    :param valuesList: lists of lists where each sub-list corresponds to a different
                       variable, and has as many values in it as the number of
                       patients (vertices)

    :param patientIdList: list of patientId values

    :return:
    """

    # PrintStartTime('Starting Multivariate PNG')

    # Create network graph with patients as the vertices and calculate edge weights based on notions of similarity
    g = graph_tool.Graph(directed=False)

    #
    # Make vertices for graph
    #

    # Iterate over all patients in patientSubset and add vertex for each
    for ii in range(len(patientIdList)):
        g.add_vertex()

    # Create vertex prop for patient ID
    patientIdProp = CreateVertexProperty(g, patientIdList, 'string')
    g.vertex_properties['idProp'] = patientIdProp

    # Normalize data and expel covariates that are uniform in value
    # print
    # print 'Starting graph construction with ' + str(len(patientIdList)) + ' vertices and ' + str(len(valuesList)) + ' covariates'
    # print '-----------------------------------------------------------------'
    nmz = list()                    # list of normalized variables
    nmzValuesList = list()          # list of
    # binary list of variables used, based on whether there is diversity in values
    varsUsedList = list()
    for ii in range(len(allVarsInclType)):
        # Check for identical values
        if len( [ jj for jj,xx in enumerate(valuesList[ii]) if xx == valuesList[ii][0]] ) == len(valuesList[ii]):
            # print 'Variable ' + str(ii) + ':  WARNING! ALL VALUES OF VARIABLE ' + allVarsInclLabels[ii] + ' ARE EQUAL'
            varsUsedList.append(0)
            # sys.exit()
        else:
            nmzValuesList.append(valuesList[ii])
            if allVarsInclType[ii] == 'continuous':
                nmz.append(PerformTStatNorm(valuesList[ii]))
                # print 'Variable ' + str(ii) + ':  ' + allVarsInclLabels[ii] + ' is continuous'
            elif allVarsInclType[ii] == 'categoricalFromContinuous':
                nmz.append(PerformTStatNorm(PerformChiSquareNorm(valuesList[ii])))
                # print 'Variable ' + str(ii) + ':  ' + allVarsInclLabels[ii] + ' is categoricalFromContinuous'
            elif allVarsInclType[ii] == 'categorical':
                chisq = PerformChiSquareNorm(valuesList[ii])
                if len( [ jj for jj,xx in enumerate(chisq) if xx == chisq[0]] ) == len(valuesList[ii]):
                    nmz.append([0.0] * len(chisq))
                else:
                    nmz.append(PerformTStatNorm(chisq))                # nmz.append(PerformTStatNorm(PerformChiSquareNorm(valuesList[ii])))
                # print 'Variable ' + str(ii) + ':  ' + allVarsInclLabels[ii] + ' is categorical'
            elif allVarsInclType[ii] == 'binary':
                chisq = PerformChiSquareNorm(valuesList[ii])
                if len( [ jj for jj,xx in enumerate(chisq) if xx == chisq[0]] ) == len(valuesList[ii]):
                    nmz.append([0.0] * len(chisq))
                else:
                    nmz.append(PerformTStatNorm(chisq))
                # nmz.append(PerformTStatNorm(PerformChiSquareNorm(valuesList[ii])))
                # print 'Variable ' + str(ii) + ':  ' + allVarsInclLabels[ii] + ' is binary'
            varsUsedList.append(1)
            # print valuesList[ii]
            # print nmz[-1]
            # print
    # print '-----------------------------------------------------------------'
    # print 'Of original ' + str(len(varsUsedList)) + ' covariates, ' + str(sum(varsUsedList)) + ' were used.'

    # Create index to track which variables are actually used to make graph
    varsUsed = g.new_graph_property('object')
    varsUsed[g] = varsUsedList
    g.graph_properties['varsUsed'] = varsUsed

    # Store all values in valuesList in vertex Prop even if they all don't get
    # used to make graph
    clinVals = CreatePatientDataVertexProperty(g, PivotList(valuesList))
    g.vertex_properties['clinVals'] = clinVals

    graphNameProp = g.new_graph_property('string')
    graphNameProp[g] = 'g0'
    g.graph_properties['graphName'] = graphNameProp

    #
    # EDGE CONSTRUCTION
    #

    # Proceed to edge construction
    if sum(varsUsedList) > 0 and len(patientIdList) > 2:
        # At least 1 variables had diversity

        # Make graph property to hold number of variables based on demo, obs, and cdp vars
        numVars = g.new_graph_property('int')
        numVars[g] = len(nmzValuesList)
        g.graph_properties['numVars'] = numVars

        # Make vertex property to hold normalized clinical values
        nmzClinVals = CreatePatientDataVertexProperty(g, PivotList(nmz))
        g.vertex_properties['nmzClinVals']  = nmzClinVals

        # Make new graph property to hold emd matrices
        emdMats = g.new_graph_property('object')

        # Create edge prop for EMD between vertices using clinical value
        emdProp, ggMatsList = CreateMultivariateEmdProperty(g)
        emdMats[g] = ggMatsList
        g.edge_properties['emd'] = emdProp
        g.graph_properties['emdMats'] = emdMats

        # Create similarity prop and emdMats prop to store EMD values for each variable
        # in a list of matrices that are easy to retrieve later
        simProp = CreateMultivariateSimProp(g, 'float', ggMatsList)
        g.edge_properties['simProp'] = simProp

        # PrintEndTime('Finished Multivariate PNG')
    else:
        # NO variables had diversity; halt construction
        print 'No variables had diversity, halting graph construction'
        ggMatsList = []
        pass

    return g, ggMatsList

def CreateVertexProperty(g, valList, propType):
    """
        Method to create a vertex property holding a value in valList for
        every vertex in the graph

    :param g:
    :param valList:
    :param propType:
    :return:
    """

    newProp = g.new_vertex_property(propType)

    for ii in range(0,g.num_vertices()):
        newProp[g.vertex(ii)] = valList[ii]

    return newProp

def CreatePatientDataVertexProperty(g, valList,type='object'):
    """

        Create a vertex property that holds heterogeneous types of data in a list

    :param g:
    :param valList:
    :return:
    """

    # Create new property to hold raw clinical values
    patientDataVertexProp = g.new_vertex_property(type)

    # Iterate thru each vertex and add raw clinical values to each node
    for ii in range(len(valList)):
        # curPatDataArray = patientDataArray()
        # curPatDataArray.dataValueList = valList[ii]
        # patientDataVertexProp[g.vertex(ii)] = curPatDataArray
        patientDataVertexProp[g.vertex(ii)] = valList[ii]

    return patientDataVertexProp

def ExtractListsFromVertices(vertexProp, g):
    """
        Method to extract the lists at each vertex of a vertex property,
        vertexProp, belonging to a graph, g, and to return a list of lists,
        where each sub-list is a list of values from each vertex

    :param vertexProp:
    :param g:
    :return:
    """

    # Make a list to hold lists from each vertex
    varList = list()

    # Iterate thru each vertex and extract list associated with vertexProp
    for ii in range(g.num_vertices()):
        varList.append(vertexProp[g.vertex(ii)])

    return varList

def PerformTStatNorm(vals):
    """
        Method to perform traditional t-stat normalization where the values in
        vals are converted to a Normal(0,1) distribution by subtracting the mean
        and dividing by the STDDEV
    """

    arrayVals = np.array(vals)
    valsNorm = (arrayVals - np.mean(arrayVals))/np.std(arrayVals)

    return list(valsNorm)

def PerformChiSquareNorm(valList):
    """

        Method to perform a chisquare calculation on every individual value in a list
        versus the remaining values in the list

    :param valList:
    :return:
    """

    chisq = list()
    for ii in range(len(valList)):


        # Made change to how chisquare is calcualted byincluding in the expected
        # set the value that is in the observed set; avoids the problem of singleton
        # occurrences in valList
        # tempValList = copy.copy(valList)
        # tempValList.remove(valList[ii])

        csq, p = CalcChiSquareP2CDistance(valList[ii], valList)
        chisq.append(csq)

    return chisq

def CalcChiSquareP2CDistance(val1, val2):
    """

        Method to calculate the chi-squared distance between a patient and a
        cohort of patients using values of a categorical variable

        Here, val1 is the value of the person, and val2 is a list of values for
        a cohort of patients

    :param val1:
    :param val2:
    :return:
    """

    # Initial check:  is val1 in val2?  If not, can't be done:
    if val1 in val2:

        # Get number of patients in val2
        N2 = len(val2)

        # Get list of possible values in val1 and val2
        possVals = list(set([val1] + val2))

        # Get tally of frequencies in val2
        tally2 = Counter()
        for ii in val2:
            tally2[ii] += 1

        # Order counts by most common ccIds
        tally2 = Counter(tally2).most_common()

        keys2 = list()
        counts2 = list()
        for ii in range(len(tally2)):
            curKey, curCount = tally2[ii]
            keys2.append(curKey)
            counts2.append(curCount)

        # # Convert tally keys to concept names
        # keys2 = Counter(tally2).keys()
        # values2 = Counter(tally2).values()

        # # Get counts from histogram and normalize by total number of ccIds
        # counts2 = list()
        # ids = list()
        # for ii in range(len(keys2)):
        #     ccId, count = tally2[ii]
        #     counts2.append(count)
        #     ids.append(ccId)

        # # Find position of val1 in keys2, which is the lost of unique values in val2
        # val1Ind = [ jj for jj,xx in enumerate(values2) if xx == val1 ]

        # Make list of same length as keys2 with a 1 in the position and 0's elsewhere
        counts1 = list()
        for ii in range(len(keys2)):
            if keys2[ii] == val1:
                counts1.append(1)
            else:
                counts1.append(0)

        # Now, counts1, and counts2 are identically ordered lists of the same length
        # that can be compared using a chi-square measure
        chisq, p = chisquare(counts1, counts2)
        x=5

    else:
        print 'Cant do chi-square of categorical variable; expected variable is not found in reference set'
        print val1

    return chisq, p

def CreateMultivariateEmdProperty(g):
    """
        Method to create an EMD property for a graph using the tStat property
        that exists at every vertex

    :param g:
    :return:
    """

    # PrintStartTime('Starting Multivariate EMD Calculation')

    numVars = g.graph_properties['numVars']
    nmzClinVals = g.vertex_properties['nmzClinVals']
    emdProp = g.new_edge_property('vector<float>')

    # Create 3d array to hold emd values
    emd3D = np.ndarray((g.num_vertices(),g.num_vertices(),numVars))

    # Iterate thru each pair of vertices and calculate the abs/diff between the respective
    # nmzClinVals pairs
    for ii in range(g.num_vertices()):
        # curNmzClinVals0 = np.array(nmzClinVals[g.vertex(ii)].dataValueList)
        curNmzClinVals0 = np.array(nmzClinVals[g.vertex(ii)])
        for jj in range(g.num_vertices()):
            newEdge = g.add_edge(g.vertex(ii), g.vertex(jj))
            # curNmzClinVals1 = np.array(nmzClinVals[g.vertex(jj)].dataValueList)
            curNmzClinVals1 = np.array(nmzClinVals[g.vertex(jj)])
            emdProp[newEdge] = list(abs(np.subtract(curNmzClinVals0, curNmzClinVals1)))
            emd3D[ii,jj,:] = abs(np.subtract(curNmzClinVals0, curNmzClinVals1))

    # PrintStartTime('Finished Multivariate EMD Calculation')

    # Convert emd3D to a list of 2D arrays
    # Create empty list of 2-D arrays
    emdMats = list()
    # emdMats = listOfArrays()
    for ii in range(numVars):
        emdMats.append(list(emd3D[:,:,ii].ravel()))

    return emdProp, emdMats

def CreateMultivariateSimProp(g, propType, emdMatsList):
    """
        Method to develop a multivariate measure of graph similarity that
        is an extension of the univariate approach.  In this case, we
        construct a quasi-Gaussian (QG) multivariate model that uses
        a covariance that is based on the EMD distances derived between
        each pair of vertices for a specific variable.

    :param g:
    :param mu:
    :param propType:
    :return:
    """

    # PrintStartTime('Started calculating multivariate similarity property')

    # # Get emd values in matrix form, by variable
    # emdMats = AssembleEmdMatricesByVar(g)
    # emdMatsList = g.graph_properties['emdMats']
    matSize = int(math.sqrt(len(emdMatsList[0])))

    # Turn emdMats into matrices again
    emdMats = list()
    for ii in range(len(emdMatsList)):
        emdMats.append(np.array(emdMatsList[ii]).reshape((matSize, matSize)))

    # Get emd values for every edge
    emd = g.edge_properties['emd']
    simProp = g.new_edge_property(propType)

    # Calculate multivariate similarity by composing a quasi-Gaussian
    # equation using emd distances for all variables
    for ii in range(g.num_vertices()):
        vertI = g.vertex(ii)
        for jj in range(ii, g.num_vertices()):
            vertJ = g.vertex(jj)

            # Get current edge
            curEdge = g.edge(vertI, vertJ)

            # Get vector of EMDs for curEdge
            curEmds = emd[curEdge]
            numVars = len(curEmds)

            # Compose covariance
            cov = np.ndarray((numVars,numVars))

            # Iterate thru curEmds to fill covariance diagonals
            for kk in range(numVars):
                # Calculate diagonal of covariance
                cov[kk,kk] = math.pow((curEmds[kk] + np.mean(emdMats[kk][ii,:])+ np.mean(emdMats[kk][jj,:]))/3.0, 2)

                # Calculate off-diagonal similarity values
                for ll in range(kk+1,numVars):
                    cov[kk,ll] = (curEmds[kk] + np.mean(emdMats[kk][ii,:])+ np.mean(emdMats[kk][jj,:])) * (curEmds[ll]  + np.mean(emdMats[ll][ii,:])+ np.mean(emdMats[ll][jj,:]))/9.0
                    cov[ll,kk] = cov[kk,ll]

            # Calculate quasi-Gaussian similarity
            simProp[curEdge] = math.exp(-1 * np.dot(np.array(curEmds), np.dot(np.linalg.pinv(cov), np.reshape(np.array(curEmds), (numVars,1)))))

    # PrintEndTime('Finished calculating multivariate similarity property')

    # return simProp, emdMats
    return simProp

def PrepareCellLineNetworkData():

    #
    # IDENTIFY CELL LINES TO BE USED IN ANALYSIS BASED ON DESIRED TISSUE TYPE
    # AND AVAILABILITY OF PHARMACOLOGY DATA, AND RETRIEVE ASSOCIATED DRUG SENSITIVITY DATA
    # FOR DRUG INDEXED BY DRUGNUM
    #


    # Save each drug response to a separate file with the same metadata
    drugSensFile = list()
    for ii in range(numDrugSensitivityTypes):
        for jj in range(numDrugs):
            drugSensFile.append(savedDataFolder + expTag + drugSensitivityFile[ii][1] + '_' + drugName[jj])

    # Read in cell line pharmacology data files to get final list of cell
    # lines that is in common among all drugs

    cellLineOrig = list()
    cellLinePh = list()
    drugSensPh = list()
    numFiles = len(drugSensFile)
    for ii in range(numFiles):
        if ii == 0:
            f0 = open(drugSensFile[0],'r')
            curDrugNamePh = cPickle.load(f0)
            cellLineOrig.append(cPickle.load(f0))
            cellLinePh = cellLineOrig[0]
            # curTissuePh = cPickle.load(f0)
            drugSensPh.append(cPickle.load(f0))
            # print(str(ii) + ' :  number of cell lines is ' + str(len(cellLinePh)))
            f0.close()
        else:
            f0 = open(drugSensFile[ii],'r')
            curDrugNamePh = cPickle.load(f0)
            cellLineOrig.append(cPickle.load(f0))
            # curTissuePh = cPickle.load(f0)
            drugSensPh.append(cPickle.load(f0))
            f0.close()

            cellLinePh = sorted(list(set(cellLinePh).intersection(set(cellLineOrig[ii]))))
            # print(str(ii) + ' :  number of cell lines is ' + str(len(cellLinePh)))


    # print 'There are ' + str(len(cellLinePh)) + ' cell lines in the original pharmacology data'
    # print

    # Reopen files and get matching cell line response values for entries in cellLinePh
    # Iterate thru each drug and extract drug sens measure ; note when an NA occurs
    # and remove them afterwards to get final COMPLETE list of cell lines

    drugSensVals = list()
    for ii in range(numDrugs):
        drugSensVals.append([])
    naInds = list()

    # Iterate over each cell line
    for ii in range(len(cellLinePh)):
        # Iterate over each drug to check that there are no NAs in all drugSensPh values
        # for the current cell line
        noNas = True
        # for jj in range(len(drugSensVals)):
        for kk in range(numDrugs):
            # Find corresponding entry in cellLineOrig[kk]
            curInd = [ mm for mm,xx in enumerate(cellLineOrig[kk]) if xx == cellLinePh[ii] ][0]
            if drugSensPh[kk][curInd] == 'NA' or drugSensPh[kk+numDrugs][curInd] == 'NA':
                noNas = False

        # If no NAs, append values to drugSensVals
        if noNas == True:
            for mm in range(numDrugs):
                # print ii,mm
                curInd = [ pp for pp,xx in enumerate(cellLineOrig[mm]) if xx == cellLinePh[ii] ][0]
                drugSensVals[mm].append([float(drugSensPh[mm][curInd]), float(drugSensPh[mm + numDrugs][curInd])])
        else:
            naInds.append(ii)

    # Remove cell lines from cellLinePh that are at naInds
    for ii in sorted(naInds, reverse=True):
        cellLinePh.pop(ii)

    # print 'There are ' + str(len(cellLinePh)) + ' cell lines after NAs are removed'
    # print

    #
    # AT THIS POINT, WE HAVE A LIST OF CELL LINES (cellLinePh) AND THE ASSOCIATED
    # SENSITIVITY VALUES FOR EACH DRUG (drugSensPh), WHERE ALL THE NAs HAVE BEEN
    # REMOVED FOR EVERY DRUG; THIS IS THE LIST OF CELL LINES WITH VALID VALUES

    # REtain only those cell lines specified by tissueType
    # cellLineInfo, cosmicId, studyAbb, cancerType = ReadCellLineInfoFile_AZD9291(cellLineInfoFile[0])
    cellLineInfo, cosmicId, studyAbb, cancerType = ReadCellLineInfoFile_Mapk(cellLineInfoFile[0])
    cellLineTissPh = list()
    studyLineAbbPh = list()
    drugSensValsTiss = list()
    for ii in range(numDrugs):
        drugSensValsTiss.append([])

    for ii in range(len(cellLinePh)):
        # Check to see if current cell line is in cellLineInfo

        # 6-9-16:  Modification to allow all cell lines to pass thru regardless of whether they appear in
        # the cellLineInfo file

        if cellLineInfo.count(cellLinePh[ii]) > 0:
            matchInd = cellLineInfo.index(cellLinePh[ii])
            if studyAbb[matchInd] in tissType or tissType[0] == 'ALL':
                cellLineTissPh.append(cellLinePh[ii])
                studyLineAbbPh.append(studyAbb[matchInd])
                for jj in range(numDrugs):
                    drugSensValsTiss[jj].append(drugSensVals[jj][ii])
        else:
            if tissType[0] == 'ALL':
                cellLineTissPh.append(cellLinePh[ii])
                studyLineAbbPh.append('MISSING')
                for jj in range(numDrugs):
                    drugSensValsTiss[jj].append(drugSensVals[jj][ii])

    # print 'There are ' + str(len(cellLineTissPh)) + ' cell lines that match the desired tissue type'
    # print

    # cellLineTissPh contains all the cell lines with the appropriate tissue type
    # drugSensValsTiss contains all the drug sensitivity data for the cell lines in cellLineTissPh


    # To facilitate matching of cell lines for pharmacology and cell lines for GAs, take
    # out dashes from cell line names

    cellLineTissPhFix = list()
    for ii in range(len(cellLineTissPh)):
        cellLineTissPhFix.append(cellLineTissPh[ii].replace('-','').upper())

    # Find list of cell lines that are in all geneAlteration files and pharmacology file

    # Iterate thru each file and pull out values, mutations, and gene Name
    # cellLineAltVal = list()
    # cellLineAlt = list()
    # geneName = list()
    # for ii in range(len(geneAlterationFiles)):
    #     # Read in mutations information
    #     f0 = open(savedDataFolder + expTag + geneAlterationFiles[ii][1],'r')
    #     cellLineAltVal.append(cPickle.load(f0))
    #     cellLineAlt.append(cPickle.load(f0))
    #     geneName.append(sorted(cPickle.load(f0)))
    #     f0.close()

    cellLinesX = list()
    drugSensValsX = list()
    studyLineAbbX = list()
    for ii in range(numDrugs):
        drugSensValsX.append([])

    # Iterate thru each cell line and see if it's in the gene alteration files
    for ii in range(len(cellLineTissPhFix)):
        found = True
        # for jj in range(len(geneAlterationFiles)):
        #     if cellLineTissPhFix[ii] not in cellLineAlt[jj]:
        #         found = False
        if found == True:
            cellLinesX.append(cellLineTissPhFix[ii])
            studyLineAbbX.append(studyLineAbbPh[ii])
            for kk in range(numDrugs):
                drugSensValsX[kk].append(drugSensValsTiss[kk][ii])
        else:
            found = True

    # print 'There are ' + str(len(cellLinesX)) + ' remaining after intersection with the cell lines for which there is GA data'
    # print

    #
    #   WRITE ALL DATA TO FILE, TO BE LATER RETRIEVED TO BUILD A GRAPH
    #

    # Save final lists of data for current processing
    f1 = open(graphsFolder + expTag + curTarg + '_' + 'RawGraphData','wb')
    cPickle.dump(cellLinesX,f1)
    cPickle.dump(drugSensValsX,f1)
    cPickle.dump([],f1)
    cPickle.dump(studyLineAbbX,f1)
    f1.close()

    return cellLinesX, drugSensValsX, studyLineAbbX

def RebuildGList():
    """

        Method to rebuild gList, which is a list that contains the original
        graph and the resulting constituent graphs from the recursive
        decomposition

    :param drugNum:
    :return:
    """

    f0 = open(savedDataFolder + expTag + curTarg + '_' + 'PartitioningResults', 'r')
    iterBinPartList = cPickle.load(f0)
    gListGraphNameList = cPickle.load(f0)
    gListNumVertsList = cPickle.load(f0)
    gListVarsUsedList = cPickle.load(f0)
    gListClosure = cPickle.load(f0)
    # iterSscoreList = cPickle.load(f0)
    f0.close()

    #
    gList = list()
    cln = load_graph(graphsFolder + expTag + 'CellLineNetwork.gt')
    gList.append(cln)
    # gList.append(load_graph(graphsFolder + 'CellLineNetwork.gt'))
    for ii in range(1,len(gListGraphNameList)):
        drugName = gListGraphNameList[ii][0:len(gListGraphNameList[ii])]
        gList.append(load_graph(graphsFolder + expTag + curTarg + '_' + drugName + 'MultivariateGraph.gt'))

    return gList, gListClosure

############################## BENCHMARKING ###############################

def PerformKMeans(sensValsList, cellLineX, K):

    from sklearn.cluster import KMeans
    from matplotlib import pyplot as plt
    from matplotlib import patches as pts

    # K = 4

    kmeans = KMeans(n_clusters=K, random_state=0).fit(PivotList(sensValsList))

    ############################################################################33

    # Organize inputs (cellLineX) by subpop

    subpopIds = list()
    for ii in range(K):
        # Find cell line ids for each value of K
        subpopIds.append(sorted([ cellLineX[jj] for jj,xx in enumerate(list(kmeans.labels_)) if xx == ii   ]))

    suffix = 'KMeans_' + str(K)
    # Open up output file and write subpop data
    with open(exportsFolder + expTag + curTarg + '_' + 'ExportedSubpopData' + suffix + '.csv', 'wb') as csvfile:
        subpopWriter = csv.writer(csvfile, delimiter=',', dialect='excel')
        # Write a row for each subpopulation
        # Format is:  subpopNum, Nsubpop, subpopIds
        for ii in range(K):
            subpopWriter.writerow([str(ii+1)] + [str(len(subpopIds[ii]))] + subpopIds[ii])

    ############################################################################33

    fullClusCtr = list(kmeans.cluster_centers_)

    ic50ClusCtrs = list()
    clusSize = list()
    for ii in range(kmeans.n_clusters):
        ic50ClusCtrs.append([ list(fullClusCtr[ii])[1], list(fullClusCtr[ii])[3] ])
        clusSize.append( len([ jj for jj,xx in enumerate(list(kmeans.labels_)) if xx == ii]))

    ic50ClusCtrs = PivotList(ic50ClusCtrs)

    clusSizeFrac = [ float(xx)/float(sum(clusSize)) for ii,xx in enumerate(clusSize) ]

    # plt.figure()
    fig,ax = plt.subplots()
    plt.plot(ic50ClusCtrs[0], ic50ClusCtrs[1], 'k*')
    plt.scatter(sensValsList[1], sensValsList[3])

    for ii in range(kmeans.n_clusters):
        circle1 = plt.Circle((ic50ClusCtrs[0][ii], ic50ClusCtrs[1][ii]), radius=clusSizeFrac[ii], color='red')
        ax.add_artist(circle1)
    plt.grid()
    # plt.xlabel('SB590885')
    # plt.ylabel('CI-1040')
    # plt.savefig(figFolder + expTag + 'SB59085_CI1040_KMeans.png')
    plt.xlabel('Afatinib')
    plt.ylabel('Selumetinib-2')
    plt.savefig(figFolder + expTag + 'SB590885_CI-1040_KMeans_' + str(K) + '.png')
    plt.close()

    print kmeans.cluster_centers_


    return

################################### RECURSIVE GRAPH PARTITIONING #####################

def CompleteRecursivelyPartitionGraph(g):
    """

        Method to recursively partition a graph-based network until stopping criterion are met

        Unlike RecursivelyPartitionGraph, this method partitions all resulting
        partitions until they meet a stopping criterion

    :param g:
    :return:
    """

    # Get total number of variables used for each graph
    numVars = g.graph_properties['numVars']

    # Set up variables to be captured per iteration
    iterCtr = 0
    iterBinPartList = list()
    iterSscoreList = list()

    # gList is the list of graphs generated from each partition
    gList = list()
    gListCtr = 0
    gList.append(g)
    gListGraphNameList = list()
    gListGraphNameList.append('g0')
    gListNumVertsList = list()
    gListNumVertsList.append(g.num_vertices())
    gListVarsUsedList = list()
    gListVarsUsedList.append([1]*numVars)

    # Set up closure list to keep track of whether a specific sub-pop
    # is closed or not; if a value is 0, then no more partitioning should occur
    gListClosure = list()
    gListClosure.append(1)

    # Perform graph partitioning until resulting graphs have too few vertices or no more
    # covariates exist that have any diversity remaining
    continuePartitioning = True
    while gListCtr < len(gList) and continuePartitioning == True:

        #
        # Perform partition
        #
        g0, g1, binPart0, sscore0, varsUsed0, varsUsed1 = PartitionGraphAndSegment(gList[gListCtr], gListGraphNameList[gListCtr])


        # 2-10-17:  modify acceptance criteria for new subpops:  the new criteria uses 2 factors:
        #           1)  is the separation between the 2 new subpops greater than the threshold:  sscoreTh
        #           2)  are the sizes of the 2 new subpops both greater than the threshold:  minSubpopSize

        if sscore0 > sscoreTh and g0.num_vertices() >= minSubpopSize and g1.num_vertices() >= minSubpopSize:
        # if sscore0 > sscoreTh:
            # Current partition is acceptable; update bookkeeping on closure status, etc.

            # Close off current graph in closure list
            gListClosure[gListCtr] = -1

            # Add graphs to gList
            gList.append(g0)
            gList.append(g1)

            # Assess closure of resulting partitions for future partitioning
            if g0.num_vertices() < minNumVerts or float(sum(varsUsed0)) == 0.0:
                gListClosure.append(0)
            else:
                gListClosure.append(1)

            if g1.num_vertices() < minNumVerts or float(sum(varsUsed1)) == 0.0:
                gListClosure.append(0)
            else:
                gListClosure.append(1)

            #
            # Keep track of bookkeeping for new graph additions
            #

            # Derive position of new graphs to be added to gList
            # This should be the length of gList - 2 and lenth of gList - 1
            newGraphInds1 = len(gList) - 2
            newGraphInds2 = len(gList) - 1

            # Keep silhouette score
            iterSscoreList.append((newGraphInds1, newGraphInds2, sscore0))

            # Append binPart0
            iterBinPartList.append((newGraphInds1, newGraphInds2, binPart0))

            # Append number of patients in resulting graphs
            gListNumVertsList.append(g0.num_vertices())
            gListNumVertsList.append(g1.num_vertices())

            # Keep list of variables used in each graph
            gListVarsUsedList.append(varsUsed0)
            gListVarsUsedList.append(varsUsed1)

            # Add graph names for new graphs
            lenGraphName = len(list(gListGraphNameList))
            graphName0 = g0.graph_properties['graphName']
            graphName1 = g1.graph_properties['graphName']
            gListGraphNameList.append(graphName0)
            gListGraphNameList.append(graphName1)

            # Save graphs
            graphName0 = g0.graph_properties['graphName']
            g0.save(graphsFolder + expTag + curTarg + '_' + graphName0 + 'MultivariateGraph.gt')
            graphName1 = g1.graph_properties['graphName']
            g1.save(graphsFolder + expTag + curTarg + '_' + graphName1 + 'MultivariateGraph.gt')

        else:
            # Current partition is unacceptable, reject partition and move forward
            # Close off current graph in closure list
            gListClosure[gListCtr] = 0

        # Increment gListCtr by 2 to account for the 2 new graphs added to gList
        # gListCtr += 2

        # Determine whether to continue partitioning based on closure
        # value of graphs remaining in gList
        # if any(gListClosure[gListCtr-1:len(gListClosure)]) == 1:
        if any(gListClosure[gListCtr+1:len(gListClosure)]) == 1:
            # There is a graph that is not closed; now find out which is the first one
            # and move gListCtr to it
            openGraphInds = [ jj for jj,xx in enumerate(gListClosure) if xx == 1 ]
            gListCtr = openGraphInds[0]
            # Advance iteration counter
            iterCtr += 1
        else:
            continuePartitioning = False

    # Save all variables describing partition
    f0 = open(savedDataFolder + expTag + curTarg + '_' + 'PartitioningResults', 'wb')
    cPickle.dump(iterBinPartList, f0)
    cPickle.dump(gListGraphNameList, f0)
    cPickle.dump(gListNumVertsList, f0)
    cPickle.dump(gListVarsUsedList, f0)
    cPickle.dump(gListClosure, f0)
    f0.close()

    return gList, gListClosure, iterBinPartList, iterSscoreList

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

    # Extract studyLineAbb property
    studyLineAbbProp = g.vertex_properties['studyLineAbb']
    studyLineAbbVals = list()
    for ii in range(g.num_vertices()):
        studyLineAbbVals.append(studyLineAbbProp[g.vertex(ii)])
    studyLineAbbVals0 = [ studyLineAbbVals[ii] for ii,xx in enumerate(binFVec) if xx == 0 ]
    studyLineAbbVals1 = [ studyLineAbbVals[ii] for ii,xx in enumerate(binFVec) if xx == 1 ]

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
    g0.graph_properties['graphName'] = graphName0
    graphName1 = g1.new_graph_property('string', parentLabel + '1')
    g1.graph_properties['graphName'] = graphName1

    # Add parent name to each graph
    parentName0 = g0.new_graph_property('string')
    parentName0[g0] = parentLabel
    g0.graph_properties['parentName'] = parentName0
    parentName1 = g1.new_graph_property('string')
    parentName1[g1] = parentLabel
    g1.graph_properties['parentName'] = parentName1

    # Add studyLineAbb property to g0 and g1
    studyLineAbbProp0 = g0.new_vertex_property('string')
    for ii in range(g0.num_vertices()):
        studyLineAbbProp0[g0.vertex(ii)] = studyLineAbbVals0[ii]
    g0.vertex_properties['studyLineAbb'] = studyLineAbbProp0
    studyLineAbbProp1 = g1.new_vertex_property('string')
    for ii in range(g1.num_vertices()):
        studyLineAbbProp1[g1.vertex(ii)] = studyLineAbbVals1[ii]
    g1.vertex_properties['studyLineAbb'] = studyLineAbbProp1

    # print
    # print 'Graph ' + graphName0[g0] + ' has ' + str(len(ids0)) + ' vertices using ' + str(numVarsUsed0) + ' variables'
    # print 'Graph ' + graphName1[g1] + ' has ' + str(len(ids1)) + ' vertices using ' + str(numVarsUsed1) + ' variables'
    # print
    # print'***************************************************************************'
    # print

    return g0, g1, binFVec, sscore, varsUsed0, varsUsed1

def PerformLaplacian(x):

    # First, if x is an adjacency, then it's diagonal values must be = 1
    for ii in range(len(x)):
        x.itemset((ii,ii), 1.0)

    # Set l to -x
    l = -x

    # Calculate sum of rows of x
    gamma = np.sum(x,axis=1)

    for ii in range(len(gamma)):
        l.itemset((ii,ii), gamma[ii] - x.item((ii,ii)))

    return l

############# TURN SEGMENTATION RESULTS INTO DISTINCT SUBPOPULATIONS #####################

def ExtractSubpopulations(gList, gListClosure):
    """

        Method to extract subpopulation information from results of
        recursive graph segmentation

    :param gList:
    :param gListClosure:
    :return:
    """

    # Rebuild gList from all the partitions
    # gList, gListClosure, iterBinPartList, iterSscoreList = RebuildGList(drugNum)

    # Iterate thru gList and extract data about vars used and verts
    numVerts = list()
    varsUsed = list()
    numVarsUsed = list()
    gName = list()
    gListInd = list()
    clinVals = list()
    cellLines = list()

    for ii in range(len(gListClosure)):
        curVarsUsed = gList[ii].graph_properties['varsUsed']
        curClinVals = list()
        # Check to see whether graph is closed or could not be split,
        # in either case, it is a graph that is a final constituent
        # graph of the population
        if gListClosure[ii] == 0:
            # Current graph is a constituent graph; so keep track of it and its parameters
            gListInd.append(ii)
            numVerts.append(gList[ii].num_vertices())
            numVarsUsed.append(sum(curVarsUsed))
            gName.append(gList[ii].graph_properties['graphName'])
            varsUsed.append(gList[ii].graph_properties['varsUsed'])

            # Get clinical values out of current graph
            clinValsProp = gList[ii].vertex_properties['clinVals']
            curClinVals = ExtractListsFromVertices(clinValsProp, gList[ii])
            # for jj in range(len(curClinValsStr)):
            #     curClinVals.append(curClinValsStr[jj].dataValueList)

            # Get names of patients in current graph and keep track of this
            idProp = gList[ii].vertex_properties['idProp']
            curPatients = ExtractListsFromVertices(idProp, gList[ii])
            clinVals.append(curClinVals)
            cellLines.append(curPatients)

        else:
            # This graph is not a constituent graph
            pass


    # Write results to file
    f = open(localResultsFolder + expTag + curTarg + '_' + 'SubpopResults', 'wb')
    cPickle.dump(numVerts, f)
    cPickle.dump(numVarsUsed, f)
    cPickle.dump(gName, f)
    cPickle.dump(varsUsed, f)
    cPickle.dump(clinVals, f)
    cPickle.dump(cellLines, f)
    cPickle.dump(gListClosure, f)
    f.close()

    return numVerts, numVarsUsed, gName, varsUsed, clinVals, cellLines, gListInd

############# EXPORT SEGMENTATION RESULTS FOR ONE DRUG PAIR #####################

def ExportSubpopResults(suffix=[]):
    """

        Method to write POET results to file

    :param suffix:
    :return:
    """

    # Rebuild gList from all the partitions
    # gList, gListClosure = RebuildGList()

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

    # Load data from composite variable analysis; most importantly
    # retrieve the indexed order of most sensitive --> most resistant drug
    f = open(resultsFolder + expTag + curTarg + '_CompositeSensitivityVariables', 'r')
    cPickle.load(f)
    sortSubpopInds = cPickle.load(f)
    f.close()

    # Get list of cell lines names in each subpop
    gListInd = [ii for ii,xx in enumerate(gListClosure) if xx == 0]
    numSubpops = len(gListInd)

    # Open up output file and write subpop data
    with open(exportsFolder + expTag + curTarg + '_' + 'ExportedSubpopData' + suffix + '.csv', 'wb') as csvfile:
        subpopWriter = csv.writer(csvfile, delimiter=',', dialect='excel')
        # subpopWriter = csv.writer(opFolder + expTag + 'ExportedSubpopData', dialect='excel', delimiter=' ')

        # Write a row for each subpopulation
        # Format is:  subpopNum, Nsubpop, subpopIds
        for ii in sortSubpopInds:
        # for ii in range(numSubpops):
            if numSubpops == 1:
                curG = load_graph(graphsFolder + expTag + curTarg + '_' + 'CellLineNetwork.gt')
            else:
                curG = load_graph(graphsFolder + expTag + curTarg + '_' + gName[ii] + 'MultivariateGraph.gt')

            tissProp = curG.vertex_properties['studyLineAbb']
            tissTypes = ExtractListsFromVertices(tissProp,curG)
            subpopWriter.writerow([str(ii+1)] + [str(numVerts[ii])] + cellLines[ii] + tissTypes)

    # Load composite variables which are organized by numDrugs x numSubpops
    f = open(resultsFolder + expTag + curTarg + '_CompositeSensitivityVariables', 'r')
    compVar = PivotList(list(cPickle.load(f)))
    f.close()

    # Export composite sensitivity variable separately
    with open(exportsFolder + expTag + curTarg + '_' + 'ExportedCompositeVariables' + suffix + '.csv', 'wb') as csvfile:
        subpopWriter = csv.writer(csvfile, delimiter=',', dialect='excel')

        # Write a row for each subpopulation
        # Format is:  subpopNum, numVerts,  numDrugs, flattened compVars for current subpop
        # for ii in range(len(numVerts)):
        for ii in sortSubpopInds:
            curSubpopStr = []
            for vv in range(numDrugs):
                for ww in range(numVerts[ii]):
                    curSubpopStr = curSubpopStr + [str(compVar[ii][vv][ww])]
            subpopWriter.writerow([str(ii+1)] + [str(numVerts[ii])] + [str(numDrugs)] + curSubpopStr)

    return cellLines

def ExportAllSubpopResults(suffix=[]):
    """

        Method to write POET results to file

    :param suffix:
    :return:
    """

    import os
    import glob
    import csv

    with open(exportsFolder + expTag + curTarg + '_' + 'ExportedAllSubpopsData' + suffix + '.csv', 'wb') as csvfile:
        subpopWriter = csv.writer(csvfile, delimiter=',', dialect='excel')
        allFiles = sorted(glob.glob(graphsFolder + expTag + curTarg + '*.gt'))

        for ii in range(len(allFiles)):

            if ii == 0:
                head,tail = os.path.split(allFiles[ii])
                slashInds = [ jj for jj,xx in enumerate(tail) if xx == '_' ]
                MInd = [ jj for jj,xx in enumerate(tail) if xx == 'M']
                curGname = 'g'
            else:
                head, tail = os.path.split(allFiles[ii])
                slashInds = [jj for jj, xx in enumerate(tail) if xx == '_']
                MInd = [jj for jj, xx in enumerate(tail) if xx == 'M']
                curGname = tail[slashInds[3]+1:MInd[1]]

            curG = load_graph(allFiles[ii])
            curNumVerts = curG.num_vertices()
            tissProp = curG.vertex_properties['studyLineAbb']
            tissTypes = ExtractListsFromVertices(tissProp, curG)
            curCellLineIdProp = curG.vertex_properties['idProp']
            curCellLineIds = ExtractListsFromVertices(curCellLineIdProp, curG)
            subpopWriter.writerow([str(ii + 1)] + [str(curGname)] + [str(curNumVerts)] + curCellLineIds + tissTypes)

    # Find all graphs in graphsFolder with curTarg and get them in a list

    # # Write results to file
    # f = open(localResultsFolder + expTag + curTarg + '_' + 'SubpopResults', 'r')
    # numVerts = cPickle.load(f)
    # numVarsUsed = cPickle.load(f)
    # gName = cPickle.load(f)
    # varsUsed = cPickle.load(f)
    # clinVals = cPickle.load(f)
    # cellLines = cPickle.load(f)
    # gListClosure = cPickle.load(f)
    # f.close()

    # # Load data from composite variable analysis; most importantly
    # # retrieve the indexed order of most sensitive --> most resistant drug
    # f = open(resultsFolder + expTag + curTarg + '_CompositeSensitivityVariables', 'r')
    # cPickle.load(f)
    # sortSubpopInds = cPickle.load(f)
    # f.close()
    #
    # # Get list of cell lines names in each subpop
    # gListInd = [ii for ii,xx in enumerate(gListClosure) if xx == 0]
    # numSubpops = len(gListInd)
    #
    # # Open up output file and write subpop data
    # with open(exportsFolder + expTag + curTarg + '_' + 'ExportedSubpopData' + suffix + '.csv', 'wb') as csvfile:
    #     subpopWriter = csv.writer(csvfile, delimiter=',', dialect='excel')
    #     # subpopWriter = csv.writer(opFolder + expTag + 'ExportedSubpopData', dialect='excel', delimiter=' ')
    #
    #     # Write a row for each subpopulation
    #     # Format is:  subpopNum, Nsubpop, subpopIds
    #     for ii in sortSubpopInds:
    #     # for ii in range(numSubpops):
    #         if numSubpops == 1:
    #             curG = load_graph(graphsFolder + expTag + curTarg + '_' + 'CellLineNetwork.gt')
    #         else:
    #             curG = load_graph(graphsFolder + expTag + curTarg + '_' + gName[ii] + 'MultivariateGraph.gt')
    #
    #         tissProp = curG.vertex_properties['studyLineAbb']
    #         tissTypes = ExtractListsFromVertices(tissProp,curG)
    #         subpopWriter.writerow([str(ii+1)] + [str(numVerts[ii])] + cellLines[ii] + tissTypes)

    # # Load composite variables which are organized by numDrugs x numSubpops
    # f = open(resultsFolder + expTag + curTarg + '_CompositeSensitivityVariables', 'r')
    # compVar = PivotList(list(cPickle.load(f)))
    # f.close()
    #
    # # Export composite sensitivity variable separately
    # with open(exportsFolder + expTag + curTarg + '_' + 'ExportedCompositeVariables' + suffix + '.csv', 'wb') as csvfile:
    #     subpopWriter = csv.writer(csvfile, delimiter=',', dialect='excel')
    #
    #     # Write a row for each subpopulation
    #     # Format is:  subpopNum, numVerts,  numDrugs, flattened compVars for current subpop
    #     # for ii in range(len(numVerts)):
    #     for ii in sortSubpopInds:
    #         curSubpopStr = []
    #         for vv in range(numDrugs):
    #             for ww in range(numVerts[ii]):
    #                 curSubpopStr = curSubpopStr + [str(compVar[ii][vv][ww])]
    #         subpopWriter.writerow([str(ii+1)] + [str(numVerts[ii])] + [str(numDrugs)] + curSubpopStr)

    return curCellLineIds

############# VISUALIZE SEGMENTATION, WITH OUTOMES #####################

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
    numVerts, numVarsUsed, leafNames, varsUsed, junk7, junk8, gListClosure = RetrieveSubpopResults(
        resultsFolderTag)

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

def DrawHierarchicalTreeWithOutcomes(filename, resultsFolderTag, outcomeStr):
    """

        Method to draw hierarchical tree

    :return: 
    """

    import create_tree_image as cti
    from cPickle import load

    # # Get saved outcome rates overall and for subpops
    # f = open(resultsFolderTag + 'SubpopOutcomes', 'r')
    # subpopReadRate = load(f)
    # overallReadRate = load(f)
    # f.close()

    # Retrieve subpop parameters
    numVerts, numVarsUsed, leafNames, varsUsed, junk7, junk8, gListClosure = RetrieveSubpopResults(
        resultsFolderTag)

    outcomeRates = [ sum(numVerts) ] + numVerts


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
    cti.create_tree(numVerts, leafNames, homoVars, filename, outcomeRates, outcomeStr)

    return

######################### SUBPOPULATION ANALYTICS ###############################

def MakeIc50Subpop(clinVals, numVerts):
    """

        Method to plot the mean values for logIc50 for each subpopulation along
        with results of significance testing between the first drug and the others

    :param meanIc50:
    :param clinVals:
    :param numVerts:
    :return:
    """

    # Extract  logIc50 data for the 3 drugs, in each subpop
    meanIc50 = [ [] for ii,xx in enumerate(range(numDrugs)) ]
    stdIc50 = [ [] for ii,xx in enumerate(range(numDrugs)) ]
    for ii in range(len(clinVals)):
        for jj in range(numDrugs):
            meanIc50[jj].append(np.mean(PivotList(clinVals[ii])[2*jj + 1]))
            stdIc50[jj].append(np.std(PivotList(clinVals[ii])[2*jj + 1]))

    # Get number of subpops
    Nsubpops = len(meanIc50[0])
    colorStr = ['b','r','g', 'm', 'c', 'k', 'y', 'b', 'r', 'g', 'm']

    # Calculate tTEsts
    tTests = list()
    for ii in range(len(clinVals)):
        curTestVals = PivotList(clinVals[ii])[1]
        curTtest = list()
        for jj in range(1,numDrugs):
            curBaseVals = PivotList(clinVals[ii])[2*jj + 1]
            ###### USE TWO SIDED TEST
            # curT, curProb = ttest(curTestVals, curBaseVals)
            # if curProb < 0.01:
            #     curTtest.append('**')
            # elif curProb < 0.05:
            #     curTtest.append('*')
            # else:
            #     curTtest.append('')
            ###### USE ONE SIDED TEST
            curT, curProb = ttest(curTestVals, curBaseVals)
            if curProb/2.0 < 0.01 and curT < 0:
                curTtest.append('**')
            elif curProb/2.0 < 0.05 and curT < 0 :
                curTtest.append('*')
            else:
                curTtest.append('')
        tTests.append(curTtest)

    # Make xtr labels
    xtickStrs = list()
    for ii in range(Nsubpops):
        xtickStrs.append(str(ii+1) + ' (' + str(numVerts[ii])+')')

    # Determine spacing parameters for numDrugs
    barWidth = 0.90 / float(numDrugs)

    # Make figure for different logIc50 responses in each subpop for each drug
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    for ii in range(numDrugs):
        ax1.bar(np.arange(0.60+barWidth*ii,Nsubpops+0.5,1), meanIc50[ii], barWidth, color=colorStr[ii], label=drugName[ii])

    plt.grid()
    ax1.set_xticks(np.arange(1,Nsubpops+1))
    ax1.set_xticklabels(xtickStrs, rotation='vertical')
    ax1.set_title('Average Log(IC50) Responses for Each Sub-population \n * = p < 0.05; ** = p < 0.01')
    ax1.set_xlabel('Sub-population Number & Size \n ' + str(sum(numVerts)) + ' Used')
    ax1.set_ylabel('Log(IC50)')

    # Add significance testing results
    x1,x2,y1,y2 = plt.axis()
    plt.axis([x1,x2,y1,y2])
    plt.axis([x1,x2,y1-1.0,y2])
    x1,x2,y1,y2 = plt.axis()
    delta = 1.0 / float(numDrugs)
    for jj in range(Nsubpops):
        for kk in range(0,numDrugs-1):
            plt.text(jj+0.80,y1+kk*delta,tTests[jj][kk],color=colorStr[kk+1])

    plt.legend()
    fig.tight_layout()
    plt.savefig(localResultsFolder + expTag + curTarg + '_'  + 'MultiDrugLogIc50.png')
    plt.close()

    return

def RetrieveSubpopResults(rfTag):
    """

        Method to extract all variables from subpopresults file and return them

    :param rfTag: 
    :return: 
    """

    from cPickle import load

    # Open subpopresults file and load variables
    f = open(rfTag + 'SubpopResults', 'r')

    numVerts = cPickle.load(f)
    numVarsUsed = cPickle.load(f)
    gName = cPickle.load(f)
    varsUsed = cPickle.load(f)
    clinVals = cPickle.load(f)
    cellLines = cPickle.load(f)
    gListClosure = cPickle.load(f)
    f.close()

    return numVerts, numVarsUsed, gName, varsUsed, clinVals, cellLines, gListClosure

def EvaluateSubpopDist():
    """

        Method to evaluate the distribution of tissues from which the cell lines
        originate from, within each subpopulation.

        This method has been altered from its original incarnation to create a
        composite variable that captures the impact of different pharmacology
        parameters in terms of a key attribute.  IN this case, the composite
        variable represents the sensitivity of a cell line to a compound by
        dividing the log IC50 value by the AUC.

        The ordering that results is critical for the filenaming and export
        of data, since the subpops are sorted by the value of the composite
        variable.  This variable, named here indsSMinMeanCompVar, is saved
        as part of the composite variables file

    :return:
    """

    # Load subpop results
    f = open(localResultsFolder + expTag + curTarg + '_'  + 'SubpopResults', 'r')
    numVerts = cPickle.load(f)
    numVarsUsed = cPickle.load(f)
    gName = cPickle.load(f)
    varsUsed = cPickle.load(f)
    clinVals = cPickle.load(f)
    cellLines = cPickle.load(f)
    gListClosure = cPickle.load(f)
    f.close()

    Nsubpops = len(numVerts)

    # Open original CLN
    origG = load_graph(graphsFolder + expTag + curTarg + '_' + 'CellLineNetwork.gt')
    origTissTypeProp = origG.vertex_properties['studyLineAbb']
    origTissTypes =  ExtractListsFromVertices(origTissTypeProp, origG)

    # Get pop-level stats
    uniqTissTypes = sorted(list(set(origTissTypes)))
    uniqTissTypesCounts = list()
    for ii in range(len(uniqTissTypes)):
        uniqTissTypesCounts.append(origTissTypes.count(uniqTissTypes[ii]))

    #
    # CALCULATE COMPOSITE VARIABLE FOR EVERY DRUG
    #

    # Calculate composite variable for each subpop
    compVar = list()
    meanCompVar = list()
    stdCompVar = list()
    for ii in range(Nsubpops):
        compVar.append([])
        meanCompVar.append([])
        stdCompVar.append([])
        for jj in range(numDrugs):
            curAUC = PivotList(clinVals[ii])[jj*2]
            curIc50 = PivotList(clinVals[ii])[jj*2 + 1]
            compVar[ii].append(np.divide(curIc50,curAUC))
            meanCompVar[ii].append(np.mean(compVar[ii][jj]))
            stdCompVar[ii].append(np.std(compVar[ii][jj]))

    # We want the display vars to be organized by 1) subpop, and 2) drug, so pivot
    compVar = PivotList(compVar)
    meanCompVar = PivotList(meanCompVar)
    stdCompVar = PivotList(stdCompVar)

    # Sort values above so that the most sensitive subpops are first; if there are more
    minMeanCompVar = [ min(xx) for ii,xx in enumerate(PivotList(meanCompVar)) ]
    sMinMeanCompVar = sorted(minMeanCompVar)
    indsSMinMeanCompVar = [ minMeanCompVar.index(xx) for ii,xx in enumerate(sMinMeanCompVar) ]

    # Use indsSMinMeanCompVar as the ordering of the original subplots, based on max sensitivity

    # Determine spacing parameters for numDrugs
    barWidth = 0.90 / float(numDrugs)
    colorStr = ['b','r','g', 'm', 'c', 'k', 'y', 'b', 'r', 'g', 'm']

    # Make figure for different logIc50 responses in each subpop for each drug
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    for ii in range(numDrugs):
        curSMeanCompVar = [ meanCompVar[ii][xx] for jj,xx in enumerate(indsSMinMeanCompVar) ]
        curStdCompVar = [ stdCompVar[ii][xx] for jj,xx in enumerate(indsSMinMeanCompVar) ]
        ax1.bar(np.arange(0.55+ii*barWidth,Nsubpops+0.55), curSMeanCompVar, barWidth, yerr=curStdCompVar, color=colorStr[ii])

    ax1.set_xticks(np.arange(1,Nsubpops+1))
    ax1.set_xticklabels(np.arange(1,Nsubpops+1))
    ax1.set_xlabel('Sub-population Number')
    ax1.set_ylabel('Avg. Composite Sensitivity \n (log(IC50)/AUC)')
    plt.grid()
    plt.savefig(figFolder + expTag + curTarg + '_' + 'CompositeDrugSensitivity')
    plt.close()

    # Save composite variables
    f = open(resultsFolder + expTag + curTarg + '_CompositeSensitivityVariables', 'wb')
    cPickle.dump(compVar, f)
    cPickle.dump(indsSMinMeanCompVar, f)
    f.close()

    # Iterate thr each subpop graph and extract composition of each subpop
    subpopTissCounts = list()
    for ii in range(Nsubpops):
        if Nsubpops == 1:
            curG = load_graph(graphsFolder + expTag + curTarg + '_' + 'CellLineNetwork.gt')
        else:
            curG = load_graph(graphsFolder + expTag + curTarg + '_' + gName[ii] + 'MultivariateGraph.gt')
        tissProp = curG.vertex_properties['studyLineAbb']
        tissTypes = ExtractListsFromVertices(tissProp,curG)
        subpopTissCounts.append([])
        for jj in range(len(uniqTissTypes)):
            subpopTissCounts[ii].append(tissTypes.count(uniqTissTypes[jj]))

    # Make plots for each subpop that demonstrate the tissue distribution in each,
    # and observing the ordering of subpops based on the most sensitive to the most
    # resistant
    # for ii in indsSMinMeanCompVar:
    for ii in range(len(indsSMinMeanCompVar)):
        # Here ii is the literal index into the subpop data in the order that the subpops were generated
        # sortSubpopNum is the logical index of the subpop, based on its ordering based on the composite variable
        # the latter should be used for "naming" the subpop in the figure titles and images since it corresponds
        # to the ordering of the subpops in the images
        dataSubpopInd = indsSMinMeanCompVar[ii]
        # Open figure
        fig = plt.figure()

        # Make bar plots of tissue distributions with sorted values
        ax1 = fig.add_subplot(211)
        plt.bar(range(len(uniqTissTypesCounts)), uniqTissTypesCounts, color='b', label='Original')
        plt.bar(range(len(uniqTissTypesCounts)), subpopTissCounts[dataSubpopInd], color='r', label='Subpop '+str(ii+1))
        ax1.set_ylabel('Counts')
        ax1.set_xlabel('Tissue Type')
        ax1.set_xticks(np.arange(0.4,len(uniqTissTypesCounts)+1))
        ax1.set_xticklabels(uniqTissTypes, rotation='vertical')
        ax1.set_title('Distribution of Tissue Types in Subpopulation ' + str(ii+1) + ' and\n Original Cell Line Network (N = ' + str(sum(numVerts)) + ')')
        plt.grid()
        plt.legend(loc=0)

        # Plot boxplots for each covariate used in graph
        ax2 = fig.add_subplot(212)
        plt.boxplot(PivotList(clinVals[ii]))
        ax2.set_xticks(np.arange(1.0,len(allVarsInclLabels)+1))
        ax2.set_xticklabels(allVarsInclLabels, rotation=20, horizontalalignment='center')
        ax2.set_ylabel('Output Values \n(AUC, log(IC50)')
        plt.grid()
        ax2.set_title('Subpop ' + str(ii+1) + ' (N = ' + str(numVerts[dataSubpopInd]) + ')')

        # Format fig and save and then close
        plt.tight_layout(pad=1.0)
        plt.savefig(localResultsFolder + expTag + curTarg + '_' + 'SubpopTissDist' + str(ii+1) + '.png')
        plt.close()

    # Make figure showing IC50 value for all different subpops
    MakeIc50Subpop([clinVals[xx] for jj,xx in enumerate(indsSMinMeanCompVar)], [numVerts[xx] for jj,xx in enumerate(indsSMinMeanCompVar)])

    return

######################### MAIN METHOD FOR SEGMENTING ONE DRUG PAIR ###############################

def CellLine_Master_Method():
    """

        Method to wrap entire POET module for cell lines into a single method

    :return:
    """

    # ###########################################################################
    # Read in pharmacology data for each compound associated with current target
    # and store results for each drug to a separate file
    # ###########################################################################

    ReadPharmacologyData_AZD9291(drugSensitivityFile[0][0], drugSensitivityFile[0][1],numDrugs)
    ReadPharmacologyData_AZD9291(drugSensitivityFile[1][0], drugSensitivityFile[1][1],numDrugs)

    # ###########################################################################
    # Prepare data for consumption by POET:  this involves curation of the pharm data
    # to remove cell lines for which there are absent values ('NA') and
    # only retaining cell lines for the desired tissue type.
    # ###########################################################################

    PrepareCellLineNetworkData()

    # ###########################################################################
    # Make cell line network and recursively partition it;  note that
    # gList is the list of graphs that were constructed during segmentation
    # gListClosure is the accompanying list whose values are either -1 or 0, where
    # 0 indicates the graphs that are the final end products (i.e., sub-populations)
    # that are to be used (values of -1 point to graphs that were subsequently
    # split again and hence, replaced by the two offspring).  The number of values
    # in gListClosure = 0 are the number of final sub-populations.
    #
    # iterBinPartList lists data about the partition and iterSscoreList gives the
    # silhouette score that resulted from each successive segmentation
    # ###########################################################################

    g = MakeCellLineNetwork()
    gList, gListClosure, iterBinPartList, iterSscoreList = CompleteRecursivelyPartitionGraph(g)

    # ###########################################################################
    # Extract subpopulations info; this turns the results of the recursive
    # partitioning above into data about sub-populations that is
    # consumed in post-processing and subsequent enrichment analysis
    # ###########################################################################

    numVerts, numVarsUsed, gName, varsUsed, clinVals, cellLines, gListInd = ExtractSubpopulations(gList, gListClosure)

    # ###########################################################################
    # Create plots for each subpopulation showing drug response ranges for
    # each compound and distribution of tissue types
    # ###########################################################################

    EvaluateSubpopDist()

    # ###########################################################################
    # Export subpopulation results to csv file
    # ###########################################################################

    ExportSubpopResults(tissType[0])
    ExportAllSubpopResults(tissType[0])

    return numVerts

def main():
    """

        The main function of this Wrapper file executes POET for a sequence of
        pairs of compounds.  For each comopound pair (CP), a set of global variables
        are set to the specifics of the current iteration.  All of the CPs are
        spelled out in a separate file.

        Before any of the CPs are processed, the files containing the genetic alterations (GAs)
        are read in and stored in a python format for easy retrieval later.  The GAs provide
        the basis for the enrichment analysis at the end.

        For each CP, a set of cell lines are found for which the pharm data is complete.
        Those cell lines provide the basis for performing POET, based on the vector
        of length 4 consisting of the log(IC50) and AUC numbers for each compound.
        The segmentation parameters are set in the accompanying infofile.  The method,
        Cellline_Master_Method(), is the actual method the performs the core functions of
        POET, including visualizing and writing the segmentations result to an output
        folder.



    :return:
    """

    # Make output folder structure in opFolder
    MakeFolder(opFolder)

    ###########################################################################
    # Read in list of targets and associated drugs
    ###########################################################################

    targList, compList = ReadTargetsFile()

    # Get number of targets
    Ntargets = len(targList)

    ###########################################################################
    # Iterate through each target/CP and run POET to stratify the sub-populations
    ###########################################################################

    # Iterate thru each set of targets --> CPs
    IterateOverTargets(targList, compList, Ntargets)

main()
