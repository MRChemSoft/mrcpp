#pragma once

#include <vector>

#include "utils/mrcpp_omp.h"
#include "macros.h"

namespace mrcpp {

#define getObjectCache(T,X) \
    ObjectCache<T> &X = ObjectCache<T>::getInstance();

#ifdef HAVE_OPENMP
#define SET_CACHE_LOCK() omp_set_lock(&this->cache_lock)
#define UNSET_CACHE_LOCK() omp_unset_lock(&this->cache_lock)
#define TEST_CACHE_LOCK() omp_test_lock(&this->cache_lock)
#else
#define SET_CACHE_LOCK()
#define UNSET_CACHE_LOCK()
#define TEST_CACHE_LOCK()
#endif

template<class T>
class ObjectCache {
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
        this->highWaterMark = 0;
        this->memLoaded = 0;
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
    ObjectCache(ObjectCache<T> const &oc) { }
    ObjectCache<T> &operator=(ObjectCache<T> const &oc) { return *this; }
#ifdef HAVE_OPENMP
    omp_lock_t cache_lock;
#endif
private:
    int highWaterMark;
    int memLoaded; ///< memory occupied by loaded objects
    std::vector<T *> objs; ///< objects store
    std::vector<int> mem; ///< mem per object
};

}
