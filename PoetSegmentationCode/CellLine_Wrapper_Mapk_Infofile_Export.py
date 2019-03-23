__author__ = 'ktch802'

# This infofile provides paths to relevant files used by CellLine_Wrapper_Mapk_Export.py
#
# It is imported on the first line of CellLine_Wrapper_Mapk_Export
#
#

# Experiment tag for file naming; this string is prepended to output files
# Change the number (20) to indicate different experiments
expTag = 'Exp100_Wrapper_Mapk_'

# Set root folders from which files and subfolders can be found
dataRoot = '/media/nkeshava/HDD1/Data/CellLinePharm/'
outputRoot = '/media/nkeshava/HDD1/Output/CellLinePharm/'

# Output folder for saving results
opFolder = outputRoot + expTag[0:-1] + '/'
figFolder = opFolder + 'Figs/'
savedDataFolder = opFolder + 'SavedData/'
graphsFolder = opFolder + 'Graphs/'
resultsFolder = opFolder + 'Results/'
exportsFolder = opFolder + 'Exports/'

# wrapperFile is a list, whose entries are the locations of csv files containing the drug pairs to be evaluated
# Successive entries in wrapperFile would result in multiple lists of drug pairs being evaluated
# We consider a single list of drug pairs
wrapperFile = list()
wrapperFile.append(dataRoot + 'Exp19_BRAF_MEK_drug_pairs.csv')

# drugSensitivityFile is a list of the locations of different files, each containing a different
# pharmacological attribute for every drug;  they can be thought of as the "features" for each drug
# drugSensitivityType records the label for each drug sensitivity feature
drugSensitivityFile = list()
drugSensitivityFile.append([dataRoot + 'public_pharmacology_AUC_alpha.csv', 'pharmDataAUC'])
drugSensitivityFile.append([dataRoot + 'public_pharmacology_IC50_alpha.csv', 'pharmDataIC50'])
drugSensitivityType = ['AUC', 'IC50']
numDrugSensitivityTypes = len(drugSensitivityType)

# cellLineInfoFile is a list whose entries are the paths to .csv files containing detailed info on each
# cell line, including its tissue type, cosmic id, short name, and cancer type
cellLineInfoFile = list()
cellLineInfoFile.append(dataRoot + 'cell_line_info_alpha.csv')

# tissType is list whose entry can be set to filter only certain tissue types to be considered
# in terms of the cell lines used
tissType = ['ALL']

# Set thresholds to constrain segmentation of populations
#
# minNumVerts:  sets the minimum number of members a subpopulation must hvae in order to be considered
# for segmentation; if the threshold is not met, then subpopulation will NOT occur on that subpopulation and
# the algorithm will move on to the next candidate subpopulation to be segmented
#
# sscoreTh:  sets the silhouette metric threshold that must be exceeded in order for a segmentation of a
# subpopulation to be retained;  silhouette metric measures the similarity of two classes (ie subpopulations);
# if the similarity is too high (or the distance is too low), the segmentation is rejected and the original parent
# subpopulation is retained
#
# minSubpopSize:  sets the minimum size of a subpopulation; if a segmentation is performed and at least one of the
# resulting subpopulations is below this value, then the proposed segmentation is rejected and the parent subpopulation
# is retained

minNumVerts = 40
sscoreTh = 0.05
minSubpopSize = 20

