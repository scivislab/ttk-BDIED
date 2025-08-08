# state file generated using paraview version 5.10.1

# uncomment the following three lines to ensure this script works in future versions
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'TTK CinemaReader'
tTKCinemaReader1 = TTKCinemaReader(registrationName='TTKCinemaReader1', DatabasePath='vertical_swap_four_clusters.cdb')

# create a new 'TTK CinemaProductReader'
tTKCinemaProductReader1 = TTKCinemaProductReader(registrationName='TTKCinemaProductReader1', Input=tTKCinemaReader1)

# create a new 'TTK MergeTree'
tTKMergeandContourTreeFTM1 = TTKMergeTree(registrationName='TTKMergeandContourTreeFTM1', Input=tTKCinemaProductReader1)
tTKMergeandContourTreeFTM1.ScalarField = ['POINTS', 'scalar']
tTKMergeandContourTreeFTM1.InputOffsetField = ['POINTS', 'scalar']
tTKMergeandContourTreeFTM1.TreeType = 'Split Tree'
tTKMergeandContourTreeFTM1.UseAllCores = 0

# find source
tTKMergeandContourTreeFTM1_1 = FindSource('TTKMergeandContourTreeFTM1')

# find source
tTKMergeandContourTreeFTM1_2 = FindSource('TTKMergeandContourTreeFTM1')

# create a new 'TTK BlockAggregator'
tTKBlockAggregator1 = TTKBlockAggregator(registrationName='TTKBlockAggregator1', Input=[tTKMergeandContourTreeFTM1, OutputPort(tTKMergeandContourTreeFTM1_1,1), OutputPort(tTKMergeandContourTreeFTM1_2,2)])
tTKBlockAggregator1.FlattenInput = 0

# create a new 'TTK MergeTreeClustering'
tTKMergeTreeClustering1 = TTKMergeTreeClustering(registrationName='TTKMergeTreeClustering1', Input=tTKBlockAggregator1,
    OptionalInputclustering=None)
tTKMergeTreeClustering1.ComputeBarycenter = 1
tTKMergeTreeClustering1.Deterministic = 1
tTKMergeTreeClustering1.DimensionSpacing = 0.002
tTKMergeTreeClustering1.Epsilon1 = 0.0
tTKMergeTreeClustering1.Epsilon2 = 100.0
tTKMergeTreeClustering1.Epsilon3 = 100.0
tTKMergeTreeClustering1.ImportantPairs = 31.0
tTKMergeTreeClustering1.ImportantPairsSpacing = 0.2
tTKMergeTreeClustering1.NonImportantPairsSpacing = 8.0
tTKMergeTreeClustering1.UseAllCores = 0

#============================================================================
# rendering
#============================================================================

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1503, 793]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [69.47561645507812, 78.22714233398438, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [69.47561645507812, 78.22714233398438, 3391.536147647599]
renderView1.CameraFocalPoint = [69.47561645507812, 78.22714233398438, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 725.4497084006886
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1503, 793)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from tTKMergeTreeClustering1
tTKMergeTreeClustering1Display = Show(tTKMergeTreeClustering1, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'Scalar'
scalarLUT = GetColorTransferFunction('Scalar')
scalarLUT.RGBPoints = [0.028, 0.231373, 0.298039, 0.752941, 0.52, 0.865003, 0.865003, 0.865003, 1.027, 0.705882, 0.0156863, 0.14902]
scalarLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Scalar'
scalarPWF = GetOpacityTransferFunction('Scalar')
scalarPWF.Points = [0.028, 0.0, 0.5, 0.0, 1.02, 1.0, 0.5, 0.0]
scalarPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
tTKMergeTreeClustering1Display.Representation = 'Surface'
tTKMergeTreeClustering1Display.ColorArrayName = ['POINTS', 'Scalar']
tTKMergeTreeClustering1Display.LookupTable = scalarLUT
tTKMergeTreeClustering1Display.LineWidth = 4.0
tTKMergeTreeClustering1Display.SelectTCoordArray = 'None'
tTKMergeTreeClustering1Display.SelectNormalArray = 'None'
tTKMergeTreeClustering1Display.SelectTangentArray = 'None'
tTKMergeTreeClustering1Display.OSPRayScaleArray = 'BranchNodeID'
tTKMergeTreeClustering1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKMergeTreeClustering1Display.SelectOrientationVectors = 'None'
tTKMergeTreeClustering1Display.ScaleFactor = 359.8905395507813
tTKMergeTreeClustering1Display.SelectScaleArray = 'BranchNodeID'
tTKMergeTreeClustering1Display.GlyphType = 'Arrow'
tTKMergeTreeClustering1Display.GlyphTableIndexArray = 'BranchNodeID'
tTKMergeTreeClustering1Display.GaussianRadius = 17.994526977539063
tTKMergeTreeClustering1Display.SetScaleArray = ['POINTS', 'BranchNodeID']
tTKMergeTreeClustering1Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKMergeTreeClustering1Display.OpacityArray = ['POINTS', 'BranchNodeID']
tTKMergeTreeClustering1Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKMergeTreeClustering1Display.DataAxesGrid = 'GridAxesRepresentation'
tTKMergeTreeClustering1Display.PolarAxes = 'PolarAxesRepresentation'
tTKMergeTreeClustering1Display.ScalarOpacityFunction = scalarPWF
tTKMergeTreeClustering1Display.ScalarOpacityUnitDistance = 941.4809873412064
tTKMergeTreeClustering1Display.OpacityArrayName = ['POINTS', 'BranchNodeID']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tTKMergeTreeClustering1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tTKMergeTreeClustering1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]

# find source
tTKMergeTreeClustering1_1 = FindSource('TTKMergeTreeClustering1')

# show data from tTKMergeTreeClustering1_1
tTKMergeTreeClustering1_1Display = Show(OutputPort(tTKMergeTreeClustering1_1, 1), renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
tTKMergeTreeClustering1_1Display.Representation = 'Surface'
tTKMergeTreeClustering1_1Display.ColorArrayName = ['POINTS', 'Scalar']
tTKMergeTreeClustering1_1Display.LookupTable = scalarLUT
tTKMergeTreeClustering1_1Display.LineWidth = 4.0
tTKMergeTreeClustering1_1Display.SelectTCoordArray = 'None'
tTKMergeTreeClustering1_1Display.SelectNormalArray = 'None'
tTKMergeTreeClustering1_1Display.SelectTangentArray = 'None'
tTKMergeTreeClustering1_1Display.OSPRayScaleArray = 'BranchNodeID'
tTKMergeTreeClustering1_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKMergeTreeClustering1_1Display.SelectOrientationVectors = 'None'
tTKMergeTreeClustering1_1Display.ScaleFactor = 14.778521537780762
tTKMergeTreeClustering1_1Display.SelectScaleArray = 'BranchNodeID'
tTKMergeTreeClustering1_1Display.GlyphType = 'Arrow'
tTKMergeTreeClustering1_1Display.GlyphTableIndexArray = 'BranchNodeID'
tTKMergeTreeClustering1_1Display.GaussianRadius = 0.7389260768890381
tTKMergeTreeClustering1_1Display.SetScaleArray = ['POINTS', 'BranchNodeID']
tTKMergeTreeClustering1_1Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKMergeTreeClustering1_1Display.OpacityArray = ['POINTS', 'BranchNodeID']
tTKMergeTreeClustering1_1Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKMergeTreeClustering1_1Display.DataAxesGrid = 'GridAxesRepresentation'
tTKMergeTreeClustering1_1Display.PolarAxes = 'PolarAxesRepresentation'
tTKMergeTreeClustering1_1Display.ScalarOpacityFunction = scalarPWF
tTKMergeTreeClustering1_1Display.ScalarOpacityUnitDistance = 37.419744034651195
tTKMergeTreeClustering1_1Display.OpacityArrayName = ['POINTS', 'BranchNodeID']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tTKMergeTreeClustering1_1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 25.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tTKMergeTreeClustering1_1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 25.0, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for scalarLUT in view renderView1
scalarLUTColorBar = GetScalarBar(scalarLUT, renderView1)
scalarLUTColorBar.Title = 'Scalar'
scalarLUTColorBar.ComponentTitle = ''

# set color bar visibility
scalarLUTColorBar.Visibility = 1

# show color legend
tTKMergeTreeClustering1Display.SetScalarBarVisibility(renderView1, True)

# show color legend
tTKMergeTreeClustering1_1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# reset view to fit data
renderView1.ResetCamera(False)

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(1503, 793)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [69.47561645507812, 78.22714233398438, 3391.536147647599]
renderView1.CameraFocalPoint = [69.47561645507812, 78.22714233398438, 0.0]
renderView1.CameraParallelScale = 725.4497084006886

# save screenshot
SaveScreenshot('vertical-swap-four-clusters.cdb/barycenter_four_clusters_wsd.png', renderView1, ImageResolution=[3006, 1586],
    TransparentBackground=0)
    
    
tTKMergeTreeClustering1.Backend = 'Path Mapping Distance (TopoInVis 2022)'

# save screenshot
SaveScreenshot('vertical-swap-four-clusters.cdb/barycenter_four_clusters_pmd.png', renderView1, ImageResolution=[3006, 1586],
    TransparentBackground=0)


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')
