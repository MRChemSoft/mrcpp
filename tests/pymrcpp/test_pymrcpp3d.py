import numpy as np
import pymrcpp as mr

from math import isclose

min_scale = -4
max_depth = 25
order = 5
prec = 1e-3

corner = np.array([-1, -1, -1])
boxes = np.array([2, 2, 2])

world = mr.BoundingBox3D(min_scale, corner, boxes)

basis = mr.InterpolatingBasis(order)

MRA = mr.MultiResolutionAnalysis3D(world, basis, max_depth)


def phi_exact(x, y, z):
    beta = 100
    alpha = (beta/np.pi)**(3/2)

    return alpha*np.exp(-beta*(x**2 + y**2 + z**2))


def v_helm(x, y, z):
    mu = 10.0
    beta = 100.0
    alpha = (beta/np.pi)**(3/2)
    coef = -6.0*beta + 4*beta**2*x**2 +\
        4*beta**2*y**2 + 4*beta**2*z**2 - mu**2

    return (-1/(4.0*np.pi))*alpha*coef*np.exp(-beta*(x**2 + y**2 + z**2))


def v_pois(x, y, z):
    beta = 100.0
    alpha = (beta/np.pi)**(3/2)
    coef = -6.0*beta + 4*beta**2*x**2 +\
        4*beta**2*y**2 + 4*beta**2*z**2

    return (-1/(4.0*np.pi))*alpha*coef*np.exp(-beta*(x**2 + y**2 + z**2))


H = mr.HelmholtzOperator(MRA, 10.0, prec)
P = mr.PoissonOperator(MRA, prec)


phi_tree = mr.FunctionTree3D(MRA)
phi_tree_pois = mr.FunctionTree3D(MRA)
v_tree = mr.FunctionTree3D(MRA)
v_tree_pois = mr.FunctionTree3D(MRA)

add_tree = mr.FunctionTree3D(MRA)
mult_tree = mr.FunctionTree3D(MRA)


mr.project(prec, v_tree, v_helm)
mr.project(prec, v_tree_pois, v_pois)

mr.apply(prec, phi_tree, H, v_tree)
mr.apply(prec, phi_tree_pois, P, v_tree_pois)


def test_IsIntWorking():
    assert isclose(1.0, phi_tree.integrate(), rel_tol=prec)


def test_BBGetScale():
    assert world.getScale() == min_scale


def test_IBGetScalingOrder():
    assert basis.getScalingOrder() == order


def test_MRAGetOrder():
    assert MRA.getOrder() == order


def test_evalf_helm():
    assert isclose(phi_tree.evalf(0, 0, 0), phi_exact(0, 0, 0), rel_tol=prec)


def test_evalf_pelm():
    assert isclose(phi_tree_pois.evalf(0, 0, 0),
                   phi_exact(0, 0, 0), rel_tol=prec)


def test_add():
    mr.add(prec/10, add_tree, 1.0, phi_tree, -1, phi_tree_pois)
    assert isclose(add_tree.evalf(0, 0, 0), 0.0, abs_tol=prec*10)


def test_multiply():
    mr.multiply(prec, mult_tree, 1, phi_tree, phi_tree_pois)
    assert isclose(mult_tree.evalf(0, 0, 0),
                   phi_exact(0, 0, 0)**2, rel_tol=prec)
