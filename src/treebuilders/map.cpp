/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#include "map.h"
#include "MapCalculator.h"
#include "MultiplicationCalculator.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "add.h"
#include "grid.h"
#include "trees/FunctionNode.h"
#include "trees/FunctionTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"
#include <Eigen/Core>

namespace mrcpp {

template <int D>
void map(double prec,
         FunctionTree<D, double> &out,
         FunctionTree<D, double> &inp,
         FMap<double, double> fmap,
         int maxIter,
         bool absPrec) {

    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, double> builder;
    WaveletAdaptor<D, double> adaptor(prec, maxScale, absPrec);
    MapCalculator<D, double> calculator(fmap, inp);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    Timer clean_t;
    inp.deleteGenerated();
    clean_t.stop();

    print::time(10, "Time transform", trans_t);
    print::time(10, "Time cleaning", clean_t);
    print::separator(10, ' ');
}

template void map<1>(double prec, FunctionTree<1, double> &out, FunctionTree<1, double> &inp, FMap<double, double> fmap, int maxIter, bool absPrec);
template void map<2>(double prec, FunctionTree<2, double> &out, FunctionTree<2, double> &inp, FMap<double, double> fmap, int maxIter, bool absPrec);
template void map<3>(double prec, FunctionTree<3, double> &out, FunctionTree<3, double> &inp, FMap<double, double> fmap, int maxIter, bool absPrec);

} // Namespace mrcpp
