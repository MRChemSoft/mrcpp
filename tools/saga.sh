#!/bin/bash -l

if [ `hostname | grep -i login | wc -l` == 0 ]; then
   echo "This script MUST be run on the LOGIN node as it requires internet access"
   exit 1
fi

mrcpp_dir="$(pwd)"
source ${mrcpp_dir}/tools/saga.env

cd ${mrcpp_dir}
build_dir=${mrcpp_dir}/build
install_dir=${mrcpp_dir}/install

if [ -d "${build_dir}" ]; then
    echo "Build directory already exists, please remove"
    exit 1
else
    ./setup --prefix=${install_dir} --omp --mpi --cxx=mpicxx ${build_dir} && \
    cd ${build_dir} && \
    make && \
    OMP_NUM_THREADS=1 ctest --output-on-failure && \
    make install
fi

exit 0
