#include "Plot.h"
#include "RepresentableFunction.h"
#include "FunctionTree.h"
#include "MathUtils.h"

using namespace std;
using namespace Eigen;


/** Plot constructor

    Upper bound is placed BEFORE lower bound in argument list .
    Arguments:
    npts:	    points in each direction, default 1000
    *b:			upper bound, default (0, 0, ... , 0)
    *a:			lower bound, default (0, 0, ... , 0)
*/
template<int D>
Plot<D>::Plot(int npts, const double *a, const double *b)
        : nPoints(npts),
          fout(0) {
    setRange(a, b);
    setSuffix(Plot<D>::Line, ".line");
    setSuffix(Plot<D>::Surface, ".surf");
    setSuffix(Plot<D>::Cube, ".cube");
    setSuffix(Plot<D>::Grid, ".grid");
}

/** Set both bounds in one go

    If either bound is a NULL pointer, its value is set to (0, 0, ..., 0)

    Arguments:
    *a:			lower bound, default (0, 0, ... , 0)
    *b:			upper bound, default (0, 0, ... , 0)
*/
template<int D>
void Plot<D>::setRange(const double *a, const double *b) {
    for (int d = 0; d < D; d++) {
        if (a == 0) {
            A[d] = 0.0;
        } else {
            A[d] = a[d];
        }
        if (b == 0) {
            B[d] = 0.0;
        } else {
            B[d] = b[d];
        }
    }
}

/** Set number of plotting points

    The number of points is restricted to be the same in all directions
*/
template<int D>
void Plot<D>::setNPoints(int npts) {
    if (npts <= 0) MSG_ERROR("Invalid number of points");
    this->nPoints = npts;
}

/** Set file extension for output file

    The file name you decide for the output will get a predefined suffix
    that differentiates between different types of plot.

    Default values:
    line:    ".line"
    surface: ".surf"
    cube:    ".cube"
    grid:    ".grid"
*/
template<int D>
void Plot<D>::setSuffix(int t, const string &s) {
    this->suffix.insert(pair<int, string>(t, s));
}

/** Parametric plot of a function

    Plots the function func parametrically between the endpoints A and B,
    to a file named fname + file extension (".line" as default). Endpoints and
    the number of points must have been set at this point (default values are
    1000 points between the boundaries of the function).
*/
template<int D>
void Plot<D>::linePlot(const RepresentableFunction<D> &func, const string &fname) {
    println(20, "----------Line Plot-----------");
    if (not verifyRange()) MSG_ERROR("Zero range");
    calcLineCoordinates();
    evaluateFunction(func);
    stringstream file;
    file << fname << this->suffix[Plot<D>::Line];
    openPlot(file.str());
    writeLineData();
    closePlot();
    printout(20, endl);
}

/** Surface plot of a function

    Plots the function func in 2D between the endpoints A and B, to a file
    named fname + file extension (".surf" as default). Endpoints and
    the number of points must have been set at this point (default values are
    1000 points between the boundaries of the function). The nPoints points
    will be distributed evenly in the two dimensions.
*/
template<int D>
void Plot<D>::surfPlot(const RepresentableFunction<D> &func, const string &fname) {
    println(20, "--------Surface Plot----------");
    if (not verifyRange()) MSG_ERROR("Zero range");
    calcSurfCoordinates();
    evaluateFunction(func);
    stringstream file;
    file << fname << this->suffix[Plot<D>::Surface];
    openPlot(file.str());
    writeSurfData();
    closePlot();
    printout(20, endl);
}

/** Cubic plot of a function

    Plots the function func in 3D between the endpoints A and B, to a file
    named fname + file extension (".cube" as default). Endpoints and
    the number of points must have been set at this point (default values are
    1000 points between the boundaries of the function). The nPoints points
    will be distributed evenly in the three dimensions.
*/
template<int D>
void Plot<D>::cubePlot(const RepresentableFunction<D> &func, const string &fname) {
    println(20, "----------Cube Plot-----------");
    if (not verifyRange()) MSG_ERROR("Zero range");
    calcCubeCoordinates();
    evaluateFunction(func);
    stringstream file;
    file << fname << this->suffix[Plot<D>::Cube];
    openPlot(file.str());
    writeCubeData();
    closePlot();
    printout(20, endl);
}

/** Parametric plot of a FunctionTree

    Plots the function func parametrically between the endpoints A and B,
    to a file named fname + file extension (".line" as default). Endpoints and
    the number of points must have been set at this point (default values are
    1000 points between the boundaries of the function).
*/
template<int D>
void Plot<D>::linePlot(FunctionTree<D> &tree, const string &fname) {
    println(20, "----------Line Plot-----------");
    if (not verifyRange()) MSG_ERROR("Zero range");
    calcLineCoordinates();
    evaluateFunction(tree);
    stringstream file;
    file << fname << this->suffix[Plot<D>::Line];
    openPlot(file.str());
    writeLineData();
    closePlot();
    printout(20, endl);
}

/** Surface plot of a function

    Plots the function func in 2D between the endpoints A and B, to a file
    named fname + file extension (".surf" as default). Endpoints and
    the number of points must have been set at this point (default values are
    1000 points between the boundaries of the function). The nPoints points
    will be distributed evenly in the two dimensions.
*/
template<int D>
void Plot<D>::surfPlot(FunctionTree<D> &tree, const string &fname) {
    println(20, "--------Surface Plot----------");
    if (not verifyRange()) MSG_ERROR("Zero range");
    calcSurfCoordinates();
    evaluateFunction(tree);
    stringstream file;
    file << fname << this->suffix[Plot<D>::Surface];
    openPlot(file.str());
    writeSurfData();
    closePlot();
    printout(20, endl);
}

/** Cubic plot of a function

    Plots the function func in 3D between the endpoints A and B, to a file
    named fname + file extension (".cube" as default). Endpoints and
    the number of points must have been set at this point (default values are
    1000 points between the boundaries of the function). The nPoints points
    will be distributed evenly in the three dimensions.
*/
template<int D>
void Plot<D>::cubePlot(FunctionTree<D> &tree, const string &fname) {
    println(20, "----------Cube Plot-----------");
    if (not verifyRange()) MSG_ERROR("Zero range");
    calcCubeCoordinates();
    evaluateFunction(tree);
    stringstream file;
    file << fname << this->suffix[Plot<D>::Cube];
    openPlot(file.str());
    writeCubeData();
    closePlot();
    printout(20, endl);
}
/** Grid plot of a MWTree

    Writes a file named fname + file extension (".grid" as default) to be read
    by geomview to visualize the grid (of endNodes)	where the multiresolution
    function is defined. In MPI, each process will write a separate file, and
    will print only	nodes owned by itself (pluss the rootNodes).
*/
template<int D>
void Plot<D>::gridPlot(const MWTree<D> &tree, const string &fname) {
    println(20, "----------Grid Plot-----------");
    stringstream file;
    file << fname << this->suffix[Plot<D>::Grid] << "." << tree.getRankId();
    openPlot(file.str());
    writeGrid(tree);
    closePlot();
    printout(20, endl);
}

/** Parametric plot of a function

    Plots the function func parametrically between the endpoints A and B,
    and returns a vector of function values. Endpoints and
    the number of points must have been set at this point (default values are
    1000 points between the boundaries of the function).
*/
template<int D>
Eigen::VectorXd &Plot<D>::linePlot(RepresentableFunction<D> &func) {
    calcLineCoordinates();
    evaluateFunction(func);
    return this->values;
}

/** Surface plot of a function

    Plots the function func in 2D between the endpoints A and B, and returns
    a vector of function values. Endpoints and
    the number of points must have been set at this point (default values are
    1000 points between the boundaries of the function). The nPoints points
    will be distributed evenly in the two dimensions.
*/
template<int D>
Eigen::VectorXd &Plot<D>::surfPlot(RepresentableFunction<D> &func) {
    calcSurfCoordinates();
    evaluateFunction(func);
    return this->values;
}

/** Cubic plot of a function

    Plots the function func in 3D between the endpoints A and B, and returns
    a vector of function values. Endpoints and
    the number of points must have been set at this point (default values are
    1000 points between the boundaries of the function). The nPoints points
    will be distributed evenly in the three dimensions.
*/
template<int D>
Eigen::VectorXd &Plot<D>::cubePlot(RepresentableFunction<D> &func) {
    calcCubeCoordinates();
    evaluateFunction(func);
    return this->values;
}

/** Calculating coordinates to be evaluated

    Generating a vector of nPoints equidistant coordinates that makes up the
    straight line between A and B in D dimensions. Coordiates are stored in
    the matrix coords for later evaluation.
*/
template<int D>
void Plot<D>::calcLineCoordinates() {
    double step[D];

    if (this->nPoints <= 0) {
        MSG_ERROR("Invalid number of points for plotting");
        return;
    }
    this->coords = MatrixXd::Zero(this->nPoints, D);
    for (int d = 0; d < D; d++) {
        step[d] = (this->B[d] - this->A[d]) / (this->nPoints + 1);
    }
    for (int i = 0; i < this->nPoints; i++) {
        for (int d = 0; d < D; d++) {
            this->coords(i, d) =  this->A[d] + (i + 1) * step[d];
        }
    }
}

template<int D>
void Plot<D>::calcSurfCoordinates() {
    NOT_IMPLEMENTED_ABORT;
//    if (D != 2) {
//        MSG_ERROR("Cannot plot planes for dim != 2!");
//        return;
//    }

//    int nPerDim = (int) floor(sqrt(this->nPoints));
//    int nRealPoints = MathUtils::ipow(nPerDim, 2);
//    this->coords = MatrixXd::Zero(nRealPoints, 2);

//    double step[2];
//    for (int d = 0; d < D; d++) {
//        step[d] = (this->B[d] - this->A[d]) / (nPerDim + 1);
//    }

//    int n = 0;
//    for (int i = 1; i <= nPerDim; i++) {
//        for (int j = 1; j <= nPerDim; j++) {
//            this->coords(n, 0) = i * step[0] + this->A[0];
//            this->coords(n, 1) = j * step[1] + this->A[1];
//            n++;
//        }
//    }
}

/** Calculating coordinates to be evaluated

    Generating a vector of nPoints coordinates equally distributed between the
    dimensions, that fills the space of a cube defined by the lower- and upper
    corner (A and B). Coordiates are stored in the matrix coords for later
    evaluation.
*/
template<>
void Plot<3>::calcCubeCoordinates() {
    int nPerDim = (int) floor(cbrt(this->nPoints));
    int nRealPoints = MathUtils::ipow(nPerDim, 3);
    this->coords = MatrixXd::Zero(nRealPoints, 3);

    double step[3];
    for (int d = 0; d < 3; d++) {
        step[d] = (this->B[d] - this->A[d]) / (nPerDim - 1);
    }

    int n = 0;
    for (int i = 0; i < nPerDim; i++) {
        for (int j = 0; j < nPerDim; j++) {
            for (int k = 0; k < nPerDim; k++) {
                this->coords(n, 0) = i * step[0] + this->A[0];
                this->coords(n, 1) = j * step[1] + this->A[1];
                this->coords(n, 2) = k * step[2] + this->A[2];
                n++;
            }
        }
    }
}

template<int D>
void Plot<D>::calcCubeCoordinates() {
    NOT_IMPLEMENTED_ABORT
}

/** Evaluating a function in a set of predfined coordinates

    Given that the set of coordinates ("coords") has been calculated, this
    routine evaluates the function in these points and stores the results
    in the vector "values".
*/
template<int D>
void Plot<D>::evaluateFunction(const RepresentableFunction<D> &func) {
    int totNPoints = this->coords.rows();
    if (not (totNPoints > 0)) {
        MSG_ERROR("Coordinates not set, cannot evaluate");
        return;
    }
    double r[D];
    this->values = VectorXd::Zero(totNPoints);
    for (int i = 0; i < totNPoints; i++) {
        for (int d = 0; d < D; d++) {
            r[d] = this->coords(i, d);
        }
        this->values[i] = func.evalf(r);
    }
}

/** Evaluating a FunctionTree in a set of predfined coordinates

    Given that the set of coordinates ("coords") has been calculated, this
    routine evaluates the function in these points and stores the results
    in the vector "values".
*/
template<int D>
void Plot<D>::evaluateFunction(FunctionTree<D> &tree) {
    int totNPoints = this->coords.rows();
    if (not (totNPoints > 0)) {
        MSG_ERROR("Coordinates not set, cannot evaluate");
        return;
    }
    double r[D];
    this->values = VectorXd::Zero(totNPoints);
    for (int i = 0; i < totNPoints; i++) {
        for (int d = 0; d < D; d++) {
            r[d] = this->coords(i, d);
        }
        this->values[i] = tree.evalf(r);
    }
}

/** Writing plot data to file

    This will write the contents of the "coords" matrix along with the function
    values to the file stream fout.	File will contain on each line the point
    number (between 0 and nPoints),	coordinates 1 through D and the function
    value.
*/
template<int D>
void Plot<D>::writeLineData() {
    ostream &o = *this->fout;
    int totNPoints = this->coords.rows();
    for (int i = 0; i < totNPoints; i++) {
        o.precision(8);
        o.setf(ios::showpoint);
        for (int d = 0; d < D; d++) {
            o << this->coords(i, d) << " ";
        }
        o.precision(12);
//        o << ", ";
        o << this->values[i];
//        o << ", ";
        o << i;
        o << endl;
    }
}

template<int D>
void Plot<D>::writeSurfData() {
    NOT_IMPLEMENTED_ABORT
}


/** Writing plot data to file

    This will write a cube file (readable by jmol) of the function values
    previously calculated (the "values" vector).
*/
template<>
void Plot<3>::writeCubeData() {
    int np = 0;
    double max = -1.0e10;
    double min = 1.0e10;
    double isoval = 0.e0;
    double step[3];
    double p = 0.0;

    ofstream &o = *this->fout;

    o << "Cube file format. Generated by MRCPP.\n" << endl;

    int nPerDim = (int) floor(cbrt(this->nPoints));
    int nRealPoints = MathUtils::ipow(nPerDim, 3);

    for (int d = 0; d < 3; d++) {
        step[d] = (this->B[d] - this->A[d]) / (nPerDim - 1);
    }

    //	"%5d %12.6f %12.6f %12.6f\n"
    o.setf(ios::scientific);
    o.precision(12);
    o << 0       << " " << 0.0     << " " << 0.0     << " " << 0.0     << endl;
    o << nPerDim << " " << step[0] << " " << 0.0     << " " << 0.0     << endl;
    o << nPerDim << " " << 0.0     << " " << step[1] << " " << 0.0     << endl;
    o << nPerDim << " " << 0.0     << " " << 0.0     << " " << step[2] << endl;

    o << endl;
    for (int n = 0; n < nRealPoints; n++) {
        o << this->values[n] << " "; //12.5E
        if (n % 6 == 5)
            o << endl;
        if (this->values[n] < min)
            min = this->values[n];
        if (this->values[n] > max)
            max = this->values[n];
        p = abs(this->values[n]);
        if (p > 1.e-4 || p < 1.e+2) {
            np += 1;
            isoval += p;
        }
    }

    isoval = isoval / np;
    println(0, "Max value:" << max);
    println(0, "Min value:" << min);
    println(0, "Isovalue: " << isoval);
}

template<int D>
void Plot<D>::writeCubeData() {
    NOT_IMPLEMENTED_ABORT
}

template<>
void Plot<3>::writeNodeGrid(const MWNode<3> &node, const string &color) {
    double origin[3] = {0, 0, 0};
    double length = pow(2.0, -node.getScale());
    ostream &o = *this->fout;

    for (int d = 0; d < 3; d++) {
        origin[d] = node.getTranslation()[d] * length;
    }
    o << origin[0] << " " <<
         origin[1] << " " <<
         origin[2] << " " <<
         color <<
         origin[0] << " " <<
         origin[1] << " " <<
         origin[2] + length << " " <<
         color <<
         origin[0] << " " <<
         origin[1] + length << " " <<
         origin[2] + length << " " <<
         color <<
         origin[0] << " " <<
         origin[1] + length << " " <<
         origin[2] <<
         color <<
         endl;

    o << origin[0] << " " <<
         origin[1] << " " <<
         origin[2] << " " <<
         color <<
         origin[0] << " " <<
         origin[1] << " " <<
         origin[2] + length << " " <<
         color <<
         origin[0] + length << " " <<
         origin[1] << " " <<
         origin[2] + length << " " <<
         color <<
         origin[0] + length << " " <<
         origin[1] << " " <<
         origin[2] <<
         color <<
         endl;
    o << origin[0] << " " <<
         origin[1] << " " <<
         origin[2] << " " <<
         color <<
         origin[0] << " " <<
         origin[1] + length << " " <<
         origin[2] << " " <<
         color <<
         origin[0] + length << " " <<
         origin[1] + length << " " <<
         origin[2] << " " <<
         color <<
         origin[0] + length << " " <<
         origin[1] << " " <<
         origin[2] <<
         color <<
         endl;

    o << origin[0] + length << " " <<
         origin[1] + length << " " <<
         origin[2] + length << " " <<
         color <<
         origin[0] + length << " " <<
         origin[1] + length << " " <<
         origin[2] << " " <<
         color <<
         origin[0] + length << " " <<
         origin[1] << " " <<
         origin[2] << " " <<
         color <<
         origin[0] + length << " " <<
         origin[1] << " " <<
         origin[2] + length <<
         color <<
         endl;

    o << origin[0] + length << " " <<
         origin[1] + length << " " <<
         origin[2] + length << " " <<
         color <<
         origin[0] + length << " " <<
         origin[1] + length << " " <<
         origin[2] << " " <<
         color <<
         origin[0] << " " <<
         origin[1] + length << " " <<
         origin[2] << " " <<
         color <<
         origin[0] << " " <<
         origin[1] + length << " " <<
         origin[2] + length <<
         color <<
         endl;

    o << origin[0] + length << " " <<
         origin[1] + length << " " <<
         origin[2] + length << " " <<
         color <<
         origin[0] + length << " " <<
         origin[1] << " " <<
         origin[2] + length << " " <<
         color <<
         origin[0] << " " <<
         origin[1] << " " <<
         origin[2] + length << " " <<
         color <<
         origin[0] << " " <<
         origin[1] + length << " " <<
         origin[2] + length <<
         color <<
         endl;
}

template<int D>
void Plot<D>::writeNodeGrid(const MWNode<D> &node, const string &color) {
    NOT_IMPLEMENTED_ABORT
}

/** Writing grid data to file

    This will write a grid file (readable by geomview) of the grid (of endNodes)
    where the multiresolution function is defined.
    Currently only working in 3D.
*/
template<>
void Plot<3>::writeGrid(const MWTree<3> &tree) {
    ostream &o = *this->fout;
    o << "CQUAD" << endl;
    o.precision(6);
    string rootColor = " 1 1 1 0 ";
    string color = " 0 0 1 1 ";
    for (int i = 0; i < tree.getRootBox().size(); i++) {
        const MWNode<3> &rootNode = tree.getRootMWNode(i);
        writeNodeGrid(rootNode, rootColor);
    }
    for (int i = 0; i < tree.getNEndNodes(); i++) {
        const MWNode<3> &node = tree.getEndMWNode(i);
        writeNodeGrid(node, color);
    }
}

template<int D>
void Plot<D>::writeGrid(const MWTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT
}

/** Opening file for output

    Opens a file output stream fout for file named fname.
*/
template<int D>
void Plot<D>::openPlot(const string &fname) {
    if (fname.empty()) {
        if (this->fout == 0) {
            MSG_ERROR("Plot file not set!");
            return;
        } else if (this->fout->fail()) {
            MSG_ERROR("Plot file not set!");
            return;
        }
    } else {
        if (this->fout != 0) {
            this->fout->close();
        }
        this->fout = &this->fstrm;
        this->fout->open(fname.c_str());
        if (this->fout->bad()) {
            MSG_ERROR("File error");
            return;
        }
    }
}

/** Closing file

    Closes the file output stream fout.
*/
template<int D>
void Plot<D>::closePlot() {
    if (this->fout != 0) this->fout->close();
    this->fout = 0;
}

/** Checks the validity of the plotting range
*/
template<int D>
bool Plot<D>::verifyRange() {
    for (int d = 0; d < D; d++) {
        if (this->A[d] > this->B[d]) {
            return false;
        }
    }
    return true;
}

template class Plot<1>;
template class Plot<2>;
template class Plot<3>;
