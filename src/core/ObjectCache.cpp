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

#include "ObjectCache.h"
#include "CrossCorrelation.h"
#include "GaussQuadrature.h"
#include "MWFilter.h"
#include "functions/LegendrePoly.h"
#include "utils/Printer.h"

namespace mrcpp {

template <class T> ObjectCache<T> &ObjectCache<T>::getInstance() {
    static ObjectCache<T> theObjectCache;
    return theObjectCache;
}

template <class T> void ObjectCache<T>::clear() {
    for (unsigned int i = 0; i < this->objs.size(); i++) {
        if (this->objs[i] != nullptr) { unload(i); }
    }
}

template <class T> void ObjectCache<T>::load(int id) {
    MSG_INFO("This routine does nothing in this class.");
}

template <class T> void ObjectCache<T>::load(int id, T *new_o, int memory) {
    if (id >= this->highWaterMark) {
        for (int i = 0; i < id - this->highWaterMark + 1; i++) {
            this->objs.push_back(nullptr);
            this->mem.push_back(0);
        }
        this->highWaterMark = id;
    }
    if (this->objs[id] != nullptr) { return; }
    this->mem[id] = memory;
    this->memLoaded += memory;
    this->objs[id] = new_o;
}

template <class T> void ObjectCache<T>::unload(int id) {
    if (id < 0 or id > this->highWaterMark) { MSG_ERROR("Id out of bounds:" << id); }
    if (this->objs[id] == nullptr) {
        MSG_WARN("Object not loaded.");
        return;
    }
    this->memLoaded -= this->mem[id];
    this->mem[id] = 0;
    delete this->objs[id];
    this->objs[id] = nullptr;
}

template <class T> T &ObjectCache<T>::get(int id) {
    if (id < 0) { MSG_ERROR("Id out of bounds:" << id); }
    if (this->objs[id] == nullptr) { MSG_ERROR("Object not loaded!"); }
    return *(this->objs[id]);
}

template <class T> bool ObjectCache<T>::hasId(int id) {
    if (id > this->highWaterMark) return false;
    if (this->objs[id] == nullptr) return false;
    return true;
}

template class ObjectCache<LegendrePoly>;
template class ObjectCache<MWFilter>;
template class ObjectCache<GaussQuadrature>;
template class ObjectCache<CrossCorrelation>;

} // namespace mrcpp
