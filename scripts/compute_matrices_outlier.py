# state file generated using paraview version 5.10.1

# uncomment the following three lines to ensure this script works in future versions
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

import numpy as np
from sklearn.manifold import TSNE
from sklearn.manifold import MDS
import matplotlib.pyplot as plt

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'TTK CinemaReader'
tTKCinemaReader1 = TTKCinemaReader(registrationName='TTKCinemaReader1', DatabasePath='vertical-swap-outlier.cdb')

# create a new 'TTK CinemaProductReader'
tTKCinemaProductReader1 = TTKCinemaProductReader(registrationName='TTKCinemaProductReader1', Input=tTKCinemaReader1)

# create a new 'TTK MergeTree'
tTKMergeTree1 = TTKMergeTree(registrationName='TTKMergeTree1', Input=tTKCinemaProductReader1)
tTKMergeTree1.ScalarField = ['POINTS', 'scalar']
tTKMergeTree1.InputOffsetField = ['POINTS', 'scalar']
tTKMergeTree1.TreeType = 'Split Tree'

# find source
tTKMergeTree1_1 = FindSource('TTKMergeTree1')

# create a new 'TTK BlockAggregator'
tTKBlockAggregator1 = TTKBlockAggregator(registrationName='TTKBlockAggregator1', Input=[tTKMergeTree1, OutputPort(tTKMergeTree1_1,1)])
tTKBlockAggregator1.FlattenInput = 0

# create a new 'TTK MergeTreeDistanceMatrix'
tTKMergeTreeDistanceMatrix1 = TTKMergeTreeDistanceMatrix(registrationName='TTKMergeTreeDistanceMatrix1', Input=tTKBlockAggregator1,
    OptionalInput=None)
tTKMergeTreeDistanceMatrix1.Backend = 'Path Mapping Distance (TopoInVis 2022)'
tTKMergeTreeDistanceMatrix1.DistanceSquareRoot = 0
tTKMergeTreeDistanceMatrix1.Epsilon1 = 0.0
tTKMergeTreeDistanceMatrix1.Epsilon2 = 100.0
tTKMergeTreeDistanceMatrix1.Epsilon3 = 100.0

# set active source
SetActiveSource(tTKMergeTreeDistanceMatrix1)

tTKMergeTreeDistanceMatrix1.PathMappingLookahead = 0

# save data
SaveData('vertical-swap.cdb/dm_pmd.csv', proxy=tTKMergeTreeDistanceMatrix1, RowDataArrays=['Tree00', 'Tree01', 'Tree02', 'Tree03', 'Tree04', 'Tree05', 'Tree06', 'Tree07', 'Tree08', 'Tree09', 'Tree10', 'Tree11', 'Tree12', 'Tree13', 'Tree14', 'Tree15', 'Tree16', 'Tree17', 'Tree18', 'Tree19', 'treeID'],
    FieldAssociation='Row Data')
    
   
matrix0 = np.loadtxt(open("vertical-swap.cdb/dm_pmd.csv", "rb"), delimiter=",", skiprows=1,usecols=range(0,20))

fig, ax = plt.subplots()
ax.imshow(matrix0, cmap='coolwarm', interpolation='nearest')
ax.set_xticks(np.arange(0,20,5))
ax.set_yticks(np.arange(0,20,5))
fig.savefig("vertical-swap.cdb/dm_pmd.pdf", transparent=True, bbox_inches="tight", pad_inches=0.1)
 


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')
