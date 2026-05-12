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

#pragma once

#include <vector>

#include "MRCPP/macros.h"
#include "utils/omp_utils.h"

namespace mrcpp {

/**
 * @def getObjectCache(T, X)
 * @brief Bind a local reference @p X to the singleton ObjectCache<T>
 *
 * @details Expands to <tt>ObjectCache<T> &X = ObjectCache<T>::getInstance()</tt>
 */
#define getObjectCache(T, X) ObjectCache<T> &X = ObjectCache<T>::getInstance();

/**
 * @class ObjectCache
 * @tparam T The object type to be cached (owned via raw pointer)
 * @brief Lightweight, integer-indexed singleton cache with optional OpenMP locking and memory accounting
 *
 * @details
 * Stores pointers to objects of type @p T in a sparse vector indexed by integer id. One global instance
 * per @p T exists via the Meyers singleton (getInstance()). Virtual hooks load(), unload(), and get() can
 * be overridden by derived caches (e.g., FilterCache, ScalingCache) for on-demand construction.
 *
 * The base class initializes an OpenMP lock for the destructor cleanup path, but the load/get/unload
 * methods themselves are not automatically locked; derived classes are responsible for guarding
 * first-time insertions with MRCPP_SET_OMP_LOCK / MRCPP_UNSET_OMP_LOCK.
 *
 * The cache takes ownership of inserted objects: unload() deletes them and clear() unloads all entries.
 * Approximate memory usage per entry is tracked in the #mem vector and summed in #memLoaded.
 * Copy and assignment are deleted to enforce singleton semantics.
 */
template <class T> class ObjectCache {
public:
    /** @brief Singleton accessor (one cache per T) */
    static ObjectCache<T> &getInstance();
    /** @brief Unload and delete all loaded objects */
    virtual void clear();

    /**
     * @brief On-demand loader hook (no-op in base class)
     * @param id Integer key for the object to load
     *
     * @details Derived caches override this to construct and insert the object at @p id
     */
    virtual void load(int id);

    /**
     * @brief Insert an already-constructed object at index @p id
     * @param id      Integer key
     * @param[in] new_o   Ownership-transferred pointer to @p T
     * @param memory  Approximate byte size (added to #memLoaded)
     *
     * @details Grows internal storage if @p id exceeds the current capacity. If an object is already
     * present at @p id the call is a no-op.
     */
    void load(int id, T *new_o, int memory);

    /**
     * @brief Remove and delete the object at @p id
     * @param id Integer key of the entry to remove
     *
     * @details Subtracts the entry's accounted memory from #memLoaded. Emits a warning if the slot
     * is already empty.
     */
    virtual void unload(int id);

    /**
     * @brief Retrieve the object at @p id by reference
     * @param id Integer key
     * @return Reference to the stored object
     *
     * @note Emits an error if @p id is negative or if the slot is empty
     */
    virtual T &get(int id);

    /**
     * @brief Test whether an object is present at @p id
     * @param id Integer key
     * @return @c true if @p id ≤ #highWaterMark and the slot is non-null
     */
    bool hasId(int id);

    /** @return Number of slots allocated (including empty/null slots) */
    int getNObjs() { return this->objs.size(); }
    /** @return Total accounted memory over loaded entries */
    int getMem() { return this->memLoaded; }
    /** @return Accounted memory for a specific @p id (0 if empty) */
    int getMem(int id) { return this->mem[id]; }

protected:
    /**
     * @brief Protected constructor; initializes slot 0, memory accounting, and the OpenMP lock
     */
    ObjectCache() {
        this->objs.push_back(nullptr);
        this->mem.push_back(0);
        MRCPP_INIT_OMP_LOCK();
    }

    /**
     * @brief Destructor; clears all entries under the OpenMP lock and destroys the lock
     */
    virtual ~ObjectCache() {
        MRCPP_SET_OMP_LOCK();
        clear();
        MRCPP_UNSET_OMP_LOCK();
        MRCPP_DESTROY_OMP_LOCK();
    }

    // Non-copyable singleton.
    ObjectCache(ObjectCache<T> const &oc) = delete;
    ObjectCache<T> &operator=(ObjectCache<T> const &oc) = delete;

#ifdef MRCPP_HAS_OMP
    /** @brief OpenMP lock for derived-class synchronized sections */
    omp_lock_t omp_lock;
#endif

private:
    int highWaterMark{0};   ///< Largest index ever stored (inclusive)
    int memLoaded{0};       ///< Total accounted memory (bytes) over all currently loaded entries
    std::vector<T *> objs;  ///< Sparse array of owned pointers; nullptr denotes an empty slot
    std::vector<int> mem;   ///< Per-slot memory estimate in bytes (0 for empty slots)
};

} // namespace mrcpp
