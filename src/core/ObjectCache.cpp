#include "ObjectCache.h"
#include "MWFilter.h"
#include "GaussQuadrature.h"
#include "CrossCorrelation.h"
#include "functions/LegendrePoly.h"
#include "utils/Printer.h"

namespace mrcpp {

template<class T>
ObjectCache<T> &ObjectCache<T>::getInstance() {
    static ObjectCache<T> theObjectCache;
    return theObjectCache;
}

template<class T>
void ObjectCache<T>::clear() {
    for (unsigned int i = 0; i < this->objs.size(); i++) {
        if (this->objs[i] != 0) {
            unload(i);
        }
    }
}

template<class T>
void ObjectCache<T>::load(int id) {
    MSG_INFO("This routine does nothing in this class.");
}

template<class T>
void ObjectCache<T>::load(int id, T *new_o, int memory) {
    if (id >= this->highWaterMark) {
        for (int i = 0; i < id - this->highWaterMark + 1; i++) {
            this->objs.push_back(0);
            this->mem.push_back(0);
        }
        this->highWaterMark = id;
    }
    if (this->objs[id] != 0) {
        return;
    }
    this->mem[id] = memory;
    this->memLoaded += memory;
    this->objs[id] = new_o;
}

template<class T>
void ObjectCache<T>::unload(int id) {
    if (id < 0 or id > this->highWaterMark) {
        MSG_ERROR("Id out of bounds:" << id);
    }
    if (this->objs[id] == 0) {
        MSG_WARN ("Object not loaded.");
        return;
    }
    this->memLoaded -= this->mem[id];
    this->mem[id] = 0;
    delete this->objs[id];
    this->objs[id] = 0;
}

template<class T>
T& ObjectCache<T>::get(int id) {
    if (id < 0) {
        MSG_ERROR("Id out of bounds:" << id);
    }
    if (this->objs[id] == 0) {
        MSG_ERROR("Object not loaded!");
    }
    return *(this->objs[id]);
}

template<class T>
bool ObjectCache<T>::hasId(int id) {
    if (id > this->highWaterMark)
        return false;
    if (this->objs[id] == 0)
        return false;
    return true;
}

template class ObjectCache<LegendrePoly>;
template class ObjectCache<MWFilter>;
template class ObjectCache<GaussQuadrature>;
template class ObjectCache<CrossCorrelation>;

} // namespace mrcpp
