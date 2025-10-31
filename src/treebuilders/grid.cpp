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

#include "grid.h"
#include "AnalyticAdaptor.h"
#include "CopyAdaptor.h"
#include "DefaultCalculator.h"
#include "SplitAdaptor.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "add.h"
#include "functions/GaussExp.h"
#include "functions/Gaussian.h"
#include "functions/function_utils.h"
#include "utils/Printer.h"

namespace mrcpp {

template <int D, typename T> void build_grid(FunctionTree<D, T> &out, int scales) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    DefaultCalculator<D, T> calculator;
    SplitAdaptor<D, T> adaptor(maxScale, true);
    for (auto n = 0; n < scales; n++) builder.build(out, calculator, adaptor, 1);
}

template <int D, typename T> void build_grid(FunctionTree<D, T> &out, const RepresentableFunction<D, T> &inp, int maxIter) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    AnalyticAdaptor<D, T> adaptor(inp, maxScale);
    DefaultCalculator<D, T> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    print::separator(10, ' ');
}

template <int D> void build_grid(FunctionTree<D> &out, const GaussExp<D> &inp, int maxIter) {
    if (!out.getMRA().getWorldBox().isPeriodic()) {
        auto maxScale = out.getMRA().getMaxScale();
        TreeBuilder<D> builder;
        DefaultCalculator<D> calculator;
        for (auto i = 0; i < inp.size(); i++) {
            AnalyticAdaptor<D> adaptor(inp.getFunc(i), maxScale);
            builder.build(out, calculator, adaptor, maxIter);
        }
    } else {
        auto period = out.getMRA().getWorldBox().getScalingFactors();
        (void)period;
        for (auto i = 0; i < inp.size(); i++) {
            auto *gauss = inp.getFunc(i).copy();
            build_grid(out, *gauss, maxIter);
            delete gauss;
        }
    }
    print::separator(10, ' ');
}

template <int D, typename T> void build_grid(FunctionTree<D, T> &out, FunctionTree<D, T> &inp, int maxIter) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    CopyAdaptor<D, T> adaptor(inp, maxScale, nullptr);
    DefaultCalculator<D, T> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    print::separator(10, ' ');
}

template <int D, typename T> void build_grid(FunctionTree<D, T> &out, FunctionTreeVector<D, T> &inp, int maxIter) {
    for (auto i = 0; i < inp.size(); i++)
        if (out.getMRA() != get_func(inp, i).getMRA()) MSG_ABORT("Incompatible MRA");

    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    CopyAdaptor<D, T> adaptor(inp, maxScale, nullptr);
    DefaultCalculator<D, T> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    print::separator(10, ' ');
}

template <int D, typename T> void build_grid(FunctionTree<D, T> &out, std::vector<FunctionTree<D, T> *> &inp, int maxIter) {
    FunctionTreeVector<D, T> inp_vec;
    for (auto *t : inp) inp_vec.push_back({1.0, t});
    build_grid(out, inp_vec, maxIter);
}

template <int D, typename T> void copy_func(FunctionTree<D, T> &out, FunctionTree<D, T> &inp) {
    FunctionTreeVector<D, T> tmp_vec;
    tmp_vec.push_back(std::make_tuple(1.0, &inp));
    add(-1.0, out, tmp_vec);
}

template <int D, typename T> void copy_grid(FunctionTree<D, T> &out, FunctionTree<D, T> &inp) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA")
    out.clear();
    build_grid(out, inp);
}

template <int D> void copy_grid(CompFunction<D> &out, CompFunction<D> &inp) {
    out.free();
    out.func_ptr->data = inp.func_ptr->data;
    out.alloc(inp.Ncomp());
    for (int i = 0; i < inp.Ncomp(); i++) {
        if (inp.isreal()) build_grid(*out.CompD[i], *inp.CompD[i]);
        if (inp.iscomplex()) build_grid(*out.CompC[i], *inp.CompC[i]);
    }
}

template <int D, typename T> void clear_grid(FunctionTree<D, T> &out) {
    TreeBuilder<D, T> builder;
    DefaultCalculator<D, T> calculator;
    builder.clear(out, calculator);
}

template <int D, typename T> int refine_grid(FunctionTree<D, T> &out, int scales) {
    auto nSplit = 0;
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    SplitAdaptor<D, T> adaptor(maxScale, true);
    for (auto n = 0; n < scales; n++) {
        nSplit += builder.split(out, adaptor, true);
    }
    return nSplit;
}

template <int D, typename T> int refine_grid(FunctionTree<D, T> &out, double prec, bool absPrec) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    WaveletAdaptor<D, T> adaptor(prec, maxScale, absPrec);
    int nSplit = builder.split(out, adaptor, true);
    return nSplit;
}

template <int D, typename T> int refine_grid(FunctionTree<D, T> &out, FunctionTree<D, T> &inp) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA")
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    CopyAdaptor<D, T> adaptor(inp, maxScale, nullptr);
    auto nSplit = builder.split(out, adaptor, true);
    return nSplit;
}

template <int D, typename T> int refine_grid(FunctionTree<D, T> &out, const RepresentableFunction<D, T> &inp) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    AnalyticAdaptor<D, T> adaptor(inp, maxScale);
    int nSplit = builder.split(out, adaptor, true);
    return nSplit;
}

template void copy_grid(CompFunction<1> &out, CompFunction<1> &inp);
template void copy_grid(CompFunction<2> &out, CompFunction<2> &inp);
template void copy_grid(CompFunction<3> &out, CompFunction<3> &inp);

template void build_grid<1, double>(FunctionTree<1, double> &out, int scales);
template void build_grid<2, double>(FunctionTree<2, double> &out, int scales);
template void build_grid<3, double>(FunctionTree<3, double> &out, int scales);
template void build_grid<1>(FunctionTree<1> &out, const GaussExp<1> &inp, int maxIter);
template void build_grid<2>(FunctionTree<2> &out, const GaussExp<2> &inp, int maxIter);
template void build_grid<3>(FunctionTree<3> &out, const GaussExp<3> &inp, int maxIter);
template void build_grid<1, double>(FunctionTree<1, double> &out, const RepresentableFunction<1, double> &inp, int maxIter);
template void build_grid<2, double>(FunctionTree<2, double> &out, const RepresentableFunction<2, double> &inp, int maxIter);
template void build_grid<3, double>(FunctionTree<3, double> &out, const RepresentableFunction<3, double> &inp, int maxIter);
template void build_grid<1, double>(FunctionTree<1, double> &out, FunctionTree<1, double> &inp, int maxIter);
template void build_grid<2, double>(FunctionTree<2, double> &out, FunctionTree<2, double> &inp, int maxIter);
template void build_grid<3, double>(FunctionTree<3, double> &out, FunctionTree<3, double> &inp, int maxIter);
template void build_grid<1, double>(FunctionTree<1, double> &out, FunctionTreeVector<1, double> &inp, int maxIter);
template void build_grid<2, double>(FunctionTree<2, double> &out, FunctionTreeVector<2, double> &inp, int maxIter);
template void build_grid<3, double>(FunctionTree<3, double> &out, FunctionTreeVector<3, double> &inp, int maxIter);
template void build_grid<1, double>(FunctionTree<1, double> &out, std::vector<FunctionTree<1, double> *> &inp, int maxIter);
template void build_grid<2, double>(FunctionTree<2, double> &out, std::vector<FunctionTree<2, double> *> &inp, int maxIter);
template void build_grid<3, double>(FunctionTree<3, double> &out, std::vector<FunctionTree<3, double> *> &inp, int maxIter);
template void copy_func<1, double>(FunctionTree<1, double> &out, FunctionTree<1, double> &inp);
template void copy_func<2, double>(FunctionTree<2, double> &out, FunctionTree<2, double> &inp);
template void copy_func<3, double>(FunctionTree<3, double> &out, FunctionTree<3, double> &inp);
template void copy_grid<1, double>(FunctionTree<1, double> &out, FunctionTree<1, double> &inp);
template void copy_grid<2, double>(FunctionTree<2, double> &out, FunctionTree<2, double> &inp);
template void copy_grid<3, double>(FunctionTree<3, double> &out, FunctionTree<3, double> &inp);
template void clear_grid<1, double>(FunctionTree<1, double> &out);
template void clear_grid<2, double>(FunctionTree<2, double> &out);
template void clear_grid<3, double>(FunctionTree<3, double> &out);
template int refine_grid<1, double>(FunctionTree<1, double> &out, int scales);
template int refine_grid<2, double>(FunctionTree<2, double> &out, int scales);
template int refine_grid<3, double>(FunctionTree<3, double> &out, int scales);
template int refine_grid<1, double>(FunctionTree<1, double> &out, double prec, bool absPrec);
template int refine_grid<2, double>(FunctionTree<2, double> &out, double prec, bool absPrec);
template int refine_grid<3, double>(FunctionTree<3, double> &out, double prec, bool absPrec);
template int refine_grid<1, double>(FunctionTree<1, double> &out, FunctionTree<1, double> &inp);
template int refine_grid<2, double>(FunctionTree<2, double> &out, FunctionTree<2, double> &inp);
template int refine_grid<3, double>(FunctionTree<3, double> &out, FunctionTree<3, double> &inp);
template int refine_grid<1, double>(FunctionTree<1, double> &out, const RepresentableFunction<1, double> &inp);
template int refine_grid<2, double>(FunctionTree<2, double> &out, const RepresentableFunction<2, double> &inp);
template int refine_grid<3, double>(FunctionTree<3, double> &out, const RepresentableFunction<3, double> &inp);

template void build_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, int scales);
template void build_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, int scales);
template void build_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, int scales);
template void build_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, const RepresentableFunction<1, ComplexDouble> &inp, int maxIter);
template void build_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, const RepresentableFunction<2, ComplexDouble> &inp, int maxIter);
template void build_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, const RepresentableFunction<3, ComplexDouble> &inp, int maxIter);
template void build_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, FunctionTree<1, ComplexDouble> &inp, int maxIter);
template void build_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, FunctionTree<2, ComplexDouble> &inp, int maxIter);
template void build_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, FunctionTree<3, ComplexDouble> &inp, int maxIter);
template void build_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, FunctionTreeVector<1, ComplexDouble> &inp, int maxIter);
template void build_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, FunctionTreeVector<2, ComplexDouble> &inp, int maxIter);
template void build_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, FunctionTreeVector<3, ComplexDouble> &inp, int maxIter);
template void build_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, std::vector<FunctionTree<1, ComplexDouble> *> &inp, int maxIter);
template void build_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, std::vector<FunctionTree<2, ComplexDouble> *> &inp, int maxIter);
template void build_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, std::vector<FunctionTree<3, ComplexDouble> *> &inp, int maxIter);
template void copy_func<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, FunctionTree<1, ComplexDouble> &inp);
template void copy_func<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, FunctionTree<2, ComplexDouble> &inp);
template void copy_func<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, FunctionTree<3, ComplexDouble> &inp);
template void copy_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, FunctionTree<1, ComplexDouble> &inp);
template void copy_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, FunctionTree<2, ComplexDouble> &inp);
template void copy_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, FunctionTree<3, ComplexDouble> &inp);
template void clear_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out);
template void clear_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out);
template void clear_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out);
template int refine_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, int scales);
template int refine_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, int scales);
template int refine_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, int scales);
template int refine_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, double prec, bool absPrec);
template int refine_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, double prec, bool absPrec);
template int refine_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, double prec, bool absPrec);
template int refine_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, FunctionTree<1, ComplexDouble> &inp);
template int refine_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, FunctionTree<2, ComplexDouble> &inp);
template int refine_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, FunctionTree<3, ComplexDouble> &inp);
template int refine_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, const RepresentableFunction<1, ComplexDouble> &inp);
template int refine_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, const RepresentableFunction<2, ComplexDouble> &inp);
template int refine_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, const RepresentableFunction<3, ComplexDouble> &inp);

} // namespace mrcpp