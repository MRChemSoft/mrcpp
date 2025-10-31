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
 * @brief Convenience macro to bind a local reference @p X to the singleton
 *        instance of ObjectCache<T>.
 *
 * Expands to:
 *   ObjectCache<T> &X = ObjectCache<T>::getInstance();
 *
 * Example:
 *   getObjectCache(MyType, cache);
 *   if (!cache.hasId(id)) cache.load(id);
 *   MyType& obj = cache.get(id);
 */
#define getObjectCache(T, X) ObjectCache<T> &X = ObjectCache<T>::getInstance();

/**
 * @class ObjectCache
 * @tparam T  The object type to be cached (owned via raw pointer).
 *
 * @brief A lightweight, index-addressed cache with singleton access,
 *        optional OpenMP locking, and simple memory accounting.
 *
 * High-level
 * ----------
 * - Stores pointers to objects of type T in a sparse, integer-indexed array.
 * - One global instance per T (Meyers singleton via getInstance()).
 * - Provides virtual hooks `load(id)`, `unload(id)`, `get(id)` for derived
 *   caches to specialize on-demand construction and retrieval.
 * - Tracks approximate memory usage per entry and in total.
 *
 * Thread-safety
 * -------------
 * - The base class initializes an OpenMP lock (if MRCPP_HAS_OMP), and its
 *   destructor clears under that lock. However, *load/get/unload* here are not
 *   automatically locked; derived classes are expected to guard first-time
 *   construction (see FilterCache, ScalingCache, etc.).
 *
 * Ownership & lifetime
 * --------------------
 * - The cache owns stored objects: `unload(id)` deletes them.
 * - `clear()` unloads all present entries.
 * - Copy/assignment are deleted to enforce singleton semantics.
 *
 * Indexing model
 * --------------
 * - `highWaterMark` tracks the largest index ever seen.
 * - Vectors `objs` and `mem` grow to accommodate new ids; gaps are filled with
 *   `nullptr` and `0`. Presence is tested with `hasId(id)`.
 *
 * Memory accounting
 * -----------------
 * - `mem[id]` holds an approximate byte size for entry `id` (provided by the
 *   caller when inserting via `load(id, T*, memory)`).
 * - `memLoaded` sums the sizes of currently loaded entries.
 */
template <class T> class ObjectCache {
public:
    /** @brief Singleton accessor (one cache per T). */
    static ObjectCache<T> &getInstance();
    /** @brief Unload and delete all loaded objects. */
    virtual void clear();

    /**
     * @brief On-demand loader hook. Default impl is a no-op; derived caches
     *        should override to construct and insert the object for @p id.
     */
    virtual void load(int id);

    /**
     * @brief Insert an already-constructed object pointer at index @p id.
     * @param id     Integer key.
     * @param new_o  Ownership-transferred pointer to T.
     * @param memory Approximate size in bytes (for accounting).
     *
     * Expands internal storage if needed. If an object is already present
     * at @p id, this is a no-op.
     */
    void load(int id, T *new_o, int memory);

    /**
     * @brief Remove and delete the object at @p id (if present).
     *        Updates memory accounting. Virtual to allow specialization.
     */
    virtual void unload(int id);

    /**
     * @brief Retrieve a reference to the loaded object at @p id.
     *        Emits errors if @p id is invalid or if no object is loaded.
     */
    virtual T &get(int id);

    /**
     * @brief Check whether an object is present at @p id.
     * @return true if id ≤ highWaterMark and objs[id] != nullptr.
     */
    bool hasId(int id);

    /** @return Number of slots allocated (including empty/null slots). */
    int getNObjs() { return this->objs.size(); }
    /** @return Total accounted memory over loaded entries. */
    int getMem() { return this->memLoaded; }
    /** @return Accounted memory for a specific @p id (0 if empty). */
    int getMem(int id) { return this->mem[id]; }

protected:
    /**
     * @brief Protected ctor initializes slot 0, memory 0, and OMP lock.
     *
     * Slot 0 is reserved/initialized so that valid ids can start at 1 if
     * desired, but the cache also happily accepts id=0.
     */
    ObjectCache() {
        this->objs.push_back(nullptr);
        this->mem.push_back(0);
        MRCPP_INIT_OMP_LOCK();
    }

    /**
     * @brief Destructor clears the cache under lock and destroys the lock.
     *
     * Ensures that concurrent threads do not race during teardown.
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
    /** @brief OpenMP lock for derived-class synchronized sections. */
    omp_lock_t omp_lock;
#endif

private:
    /** @brief Largest index ever used (inclusive). */
    int highWaterMark{0};
    /** @brief Sum of accounted memory over loaded entries. */
    int memLoaded{0};      ///< memory occupied by loaded objects
    /** @brief Sparse vector of owned pointers; nullptr denotes empty slot. */
    std::vector<T *> objs; ///< objects store
    /** @brief Per-slot memory accounting (0 if empty). */
    std::vector<int> mem;  ///< mem per object
};

} // namespace mrcpp
