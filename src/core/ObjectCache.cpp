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

/*
 * Overview
 * --------
 * This file provides the generic implementation of a very simple object cache:
 *    ObjectCache<T>
 *
 * The cache stores pointers to objects of type T, indexed by an integer id.
 * It supports:
 *   - on-demand loading (either overridden in a derived cache, or via
 *     load(id, T*, memory) with an already-constructed object),
 *   - unloading (delete + accounting),
 *   - querying if an id is present, and
 *   - retrieving a reference to a loaded object.
 *
 * Key properties:
 *   • Sparse index space:
 *       Internally, `objs` and `mem` are vectors. The `highWaterMark` records
 *       the highest id seen so far, and the vectors are expanded with nullptr/0
 *       as needed in `load(id, T*, memory)`.
 *
 *   • Memory accounting:
 *       `mem[id]` stores a byte estimate for the object at index `id`
 *       (provided by the caller). `memLoaded` accumulates the total over
 *       loaded entries. There is no automatic eviction policy here; derived
 *       caches may use these numbers for their own management.
 *
 *   • Thread-safety:
 *       This base class does not synchronize access. Derived caches (e.g.,
 *       filter caches) add OpenMP locks around load/insert to ensure safety.
 *
 *   • Lifetime:
 *       Objects are owned by the cache (deleted on unload). `clear()` iterates
 *       over all indices and unloads any present objects.
 *
 * Explicit instantiations at the end make sure the compiler emits code for
 * the most common cached types used within MRCPP.
 */

#include "ObjectCache.h"
#include "CrossCorrelation.h"
#include "GaussQuadrature.h"
#include "MWFilter.h"
#include "functions/LegendrePoly.h"
#include "utils/Printer.h"

namespace mrcpp {

/*
 * getInstance()
 * -------------
 * Meyers' singleton accessor for ObjectCache<T>.
 * A single cache instance per T exists process-wide.
 */
template <class T> ObjectCache<T> &ObjectCache<T>::getInstance() {
    static ObjectCache<T> theObjectCache;
    return theObjectCache;
}

/*
 * clear()
 * -------
 * Unload all currently loaded objects by iterating the index range and
 * calling unload(i) for each non-null entry.
 */
template <class T> void ObjectCache<T>::clear() {
    for (unsigned int i = 0; i < this->objs.size(); i++) {
        if (this->objs[i] != nullptr) { unload(i); }
    }
}

/*
 * load(id)
 * --------
 * Default "do nothing" loader. The intent is that specialized caches
 * (e.g., FilterCache, CrossCorrelationCache) override this method to
 * construct/load the appropriate object for the given id. Calling this
 * base implementation only prints an info message.
 */
template <class T> void ObjectCache<T>::load(int id) {
    MSG_INFO("This routine does nothing in this class.");
}

/*
 * load(id, new_o, memory)
 * -----------------------
 * Insert an already-constructed object pointer at index `id`.
 * - Expands internal storage if `id` exceeds the current highWaterMark,
 *   filling with nullptr/0.
 * - If an object is already present at `id`, the call is a no-op.
 * - Otherwise, records the memory estimate, updates `memLoaded`, and stores
 *   the pointer. Ownership is transferred to the cache (deleted in unload()).
 */
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

/*
 * unload(id)
 * ----------
 * Remove and delete the object at index `id`, updating memory accounting.
 * - Validates bounds.
 * - Warns (and returns) if the slot is already empty.
 * - Sets the slot to nullptr and zeroes its memory entry.
 */
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

/*
 * get(id)
 * -------
 * Return a reference to the object stored at `id`.
 * - Emits an error if `id` is negative or if the object is not loaded.
 *   (Note: derived caches typically call hasId()/load() to ensure presence.)
 */
template <class T> T &ObjectCache<T>::get(int id) {
    if (id < 0) { MSG_ERROR("Id out of bounds:" << id); }
    if (this->objs[id] == nullptr) { MSG_ERROR("Object not loaded!"); }
    return *(this->objs[id]);
}

/*
 * hasId(id)
 * ---------
 * Query whether an object for `id` is present in the cache.
 * Returns false if `id` exceeds the current high-water mark or if the
 * slot holds nullptr; true otherwise.
 */
template <class T> bool ObjectCache<T>::hasId(int id) {
    if (id > this->highWaterMark) return false;
    if (this->objs[id] == nullptr) return false;
    return true;
}

/*
 * Explicit template instantiations
 * --------------------------------
 * Force code generation for these commonly cached types, so users linking
 * to MRCPP do not need to instantiate ObjectCache<T> themselves.
 */
template class ObjectCache<LegendrePoly>;
template class ObjectCache<MWFilter>;
template class ObjectCache<GaussQuadrature>;
template class ObjectCache<CrossCorrelation>;

} // namespace mrcpp