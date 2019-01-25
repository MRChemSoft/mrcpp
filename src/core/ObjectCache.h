/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#pragma once

#include <vector>

#include "MRCPP/macros.h"
#include "utils/omp_utils.h"

namespace mrcpp {

#define getObjectCache(T, X) ObjectCache<T> &X = ObjectCache<T>::getInstance();

#ifdef HAVE_OPENMP
#define SET_CACHE_LOCK() omp_set_lock(&this->cache_lock)
#define UNSET_CACHE_LOCK() omp_unset_lock(&this->cache_lock)
#define TEST_CACHE_LOCK() omp_test_lock(&this->cache_lock)
#else
#define SET_CACHE_LOCK()
#define UNSET_CACHE_LOCK()
#define TEST_CACHE_LOCK()
#endif

template <class T> class ObjectCache {
public:
    static ObjectCache<T> &getInstance();
    virtual void clear();

    virtual void load(int id);
    void load(int id, T *new_o, int memory);
    virtual void unload(int id);

    virtual T &get(int id);
    bool hasId(int id);

    int getNObjs() { return this->objs.size(); }
    int getMem() { return this->memLoaded; }
    int getMem(int id) { return this->mem[id]; }

protected:
    ObjectCache() {
        this->objs.push_back(0);
        this->mem.push_back(0);
#ifdef HAVE_OPENMP
        omp_init_lock(&cache_lock);
#endif
    }

    virtual ~ObjectCache() {
        SET_CACHE_LOCK();
        clear();
        UNSET_CACHE_LOCK();
#ifdef HAVE_OPENMP
        omp_destroy_lock(&cache_lock);
#endif
    }
    ObjectCache(ObjectCache<T> const &oc) = delete;
    ObjectCache<T> &operator=(ObjectCache<T> const &oc) = delete;
#ifdef HAVE_OPENMP
    omp_lock_t cache_lock;
#endif
private:
    int highWaterMark{0};
    int memLoaded{0};      ///< memory occupied by loaded objects
    std::vector<T *> objs; ///< objects store
    std::vector<int> mem;  ///< mem per object
};

} // namespace mrcpp
