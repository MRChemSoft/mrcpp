import numpy as np
import pymrcpp as mr

from math import isclose

min_scale = -4
max_depth = 25
order = 7
prec = 1e-5

corner = -1
boxes = 2

world = mr.BoundingBox1D(min_scale, corner, boxes)

basis = mr.InterpolatingBasis(order)

MRA = mr.MultiResolutionAnalysis1D(world, basis, max_depth)


def phi(x):
    beta = 100
    alpha = (beta/np.pi)**(1/2)

    return alpha*np.exp(-beta*x**2)


phi_tree = mr.FunctionTree1D(MRA)
add_tree = mr.FunctionTree1D(MRA)
mult_tree = mr.FunctionTree1D(MRA)

mr.project(prec, phi_tree, phi)


def test_IsIntWorking():
    assert isclose(1.0, phi_tree.integrate(), rel_tol=prec)


def test_BBGetScale():
    assert world.getScale() == min_scale


def test_IBGetScalingOrder():
    assert basis.getScalingOrder() == order


def test_MRAGetOrder():
    assert MRA.getOrder() == order


def test_add():
    mr.add(prec/10, add_tree, 1.0, phi_tree, -1, phi_tree)
    assert isclose(add_tree.evalf(0), 0.0, abs_tol=prec*10)


def test_multiply():
    mr.multiply(prec, mult_tree, 1, phi_tree, phi_tree)
    assert isclose(mult_tree.evalf(0), phi(0)**2, rel_tol=prec)
