#    # load Intel compilers
#    set +x
#    source /theoryfs2/common/software/intel2018/bin/compilervars.sh intel64
#    set -x
# error with Intel 2018.3, 2019.4
#
#    ALLOPTS="-gnu-prefix=${HOST}- ${OPTS}"

# configure
${BUILD_PREFIX}/bin/cmake \
        -H${SRC_DIR} \
        -Bbuild \
        -G"Ninja" \
        -DCMAKE_INSTALL_PREFIX=${PREFIX} \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_COMPILER=${CXX} \
        -DCMAKE_PREFIX_PATH="${PREFIX}" \
        -DBUILD_STATIC_LIBS=False \
        -DENABLE_TESTS=True
        # -DENABLE_OPENMP=ON

# build
cd build
ninja -j${CPU_COUNT}

# test
ctest -j${CPU_COUNT} --output-on-failure --verbose

# install
ninja install
