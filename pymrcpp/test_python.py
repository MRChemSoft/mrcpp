import numpy as np
import pymrcpp as mr 

from math import isclose

min_scale = -4
max_depth = 25 
order = 7
prec = 1e-5

corner = np.array([-1, -1, -1])
boxes = np.array([2, 2 , 2])

world = mr.BoundingBox3D(min_scale, corner, boxes)

basis = mr.InterpolatingBasis(order)

MRA = mr.MultiResolutionAnalysis3D(world, basis, max_depth)

def phi_exact(x, y, z):
    beta = 100
    alpha = (beta/np.pi)**(3/2)

    return alpha*np.exp(-beta*(x**2 + y**2 + z**2))


phi_tree = mr.FunctionTree3D(MRA)
mr.project(prec, phi_tree, phi_exact)


def test_IsIntWorking():
    assert isclose(1.0, phi_tree.integrate(), rel_tol = prec)


def test_BBGetScale():
    assert world.getScale() == min_scale


def test_IBGetScalingOrder():
    assert basis.getScalingOrder() == order


def test_MRAGetOrder():
    assert MRA.getOrder() == order


def test_evalf():
    assert isclose(phi_tree.evalf(.2, .2 ,.2), phi_exact(.2, .2, .2), rel_tol = prec)
