rm -rf build
mkdir build
pushd build

if [ "${SHORT_OS_STR}" == "Darwin" ]; then
    CMAKE_OSX_DEPLOYMENT_TARGET="-DCMAKE_OSX_DEPLOYMENT_TARGET="
else
    CMAKE_OSX_DEPLOYMENT_TARGET=""
fi


cmake -G"$CMAKE_GENERATOR" ${CMAKE_OSX_DEPLOYMENT_TARGET} ../
make
mkdir -p $PREFIX/bin
cp chain $PREFIX/bin
popd
cp CurrentRep.txt $PREFIX
cp -r bill_plans $PREFIX
cp *.csv $PREFIX
cp -R scripts $PREFIX
