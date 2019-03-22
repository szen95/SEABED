import create_tree_image as cti
import numpy as np
import pandas as pd

# import data to visualize hirarchical tree after segmentation
data = 'Exp12_Wrapper_Mapk_1_ExportedAllSubpopsDataALL.csv'
df1 = pd.read_csv(data, index_col = 0, error_bad_lines=False)

# obtain leafNames and numVerts from dataset
df2 = df1['leafName']
df3 = df1['numVerts']

# get the leafNames and numVerts in a list
leafNames = df2.values
numVerts = df3.values
# print leafNames
# print numVerts
# print len(leafNames)
# print len(numVerts)

# list of subpopulation sizes
numVerts = [745, 501, 408, 258, 244, 226, 194, 160, 150, 144, 127, 124, 104,  99,  93,  84,  70,  57,
  56,  54,  51,  50,  48,  45,  45,  43, 36,  33,  32,  32,  28,  27,  24,  23,  21,  20,
  20]

# list of strings where each string has the format gXXX, where the list is for the leaves (final subpops) of the tree in the same order as numverts
leafNames = ['g', 'g0', 'g00', 'g000', 'g1', 'g0000', 'g00000', 'g10', 'g001', 'g000000',
 'g0010', 'g100', 'g1000', 'g0000000', 'g01', 'g11', 'g00100', 'g00101', 'g10000',
 'g00000000', 'g110', 'g000001', 'g10001', 'g00000001', 'g0000001', 'g001000',
 'g101', 'g111', 'g00001', 'g0001', 'g100010', 'g001001', 'g000000010', 'g0011',
 'g000000011', 'g100011', 'g1001']

# a list of 1s and 0s the same length as numverts where a 1 indicates that the corresponding subpopulation has a homogeneous variable (not that relevant for our this and you can make the list completely of 0s)
homoVars = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# print len(homoVars)

# the string for the colorbar
outcomeStr = 'Subpopulations'
filename = 'tree_diagram.pdf'

cti.create_tree(numVerts, leafNames, homoVars, filename, outcomeStr=outcomeStr)
