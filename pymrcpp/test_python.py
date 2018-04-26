import numpy as np
import pymrcpp as mr 

min_scale = -4
max_depth = 25 
order = 7
prec = 1e-5

corner = np.array([-1, -1, -1])
boxes = np.array([2, 2 , 2])

world = mr.BoundingBox3D(min_scale, corner, boxes)

basis = mr.InterpolatingBasis(order)

MRA = mr.MultiResolutionAnalysis3D(world, basis, max_depth)

def test_BBGetScale():
    assert world.getScale() == min_scale

def test_IBGetScalingOrder():
    assert basis.getScalingOrder() == order

def test_MRAGetOrder():
    assert MRA.getOrder() == order
