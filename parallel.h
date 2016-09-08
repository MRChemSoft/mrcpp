/*
 *
 *
 *  \date Sep 27, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif Convenience definitions for parallel builds
 */

#ifndef PARALLEL_H_
#define PARALLEL_H_

#include "config.h"

#ifdef HAVE_OPENMP

#define OPENMP
#include <omp.h>

#else

#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_set_dynamic(n)
#define omp_set_lock(x)
#define omp_unset_lock(x)
#define omp_test_lock(x)

#endif

#endif /* PARALLEL_H_ */
