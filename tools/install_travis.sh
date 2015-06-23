#!/bin/bash
#set -e

LOCAL=$PWD
mkdir -p $LOCAL/src

export LD_LIBRARY_PATH=$LOCAL/lib:$LOCAL/lib64:$LOCAL/usr/lib:$LOCAL/usr/lib64
export PKG_CONFIG_PATH=$LOCAL/lib/pkgconfig:$LOCAL/usr/lib/pkgconfig
export PATH=/usr/lib/ccache:$PATH:$LOCAL/bin:$LOCAL/usr/bin
export HDF5_DIR=$LOCAL/usr

ls /usr/lib/ccache

apt-get download \
libfftw3-3 \
libfftw3-dev \
libhdf5-serial-1.8.4 \
libhdf5-serial-dev \
libblas3gf \
libblas-dev \
liblapack3gf \
liblapack-dev \
libgsl0ldbl \
libgsl0-dev

for PKG in *.deb; do
    dpkg-deb -x $PKG $LOCAL
done

wget http://gfortran.com/download/x86_64/gfortran-fetch-and-install.sh
sh gfortran-fetch-and-install.sh -d $PWD 4.7

echo usr/include
pip install h5py

# install the version of swig that for some reason we have to use
wget http://downloads.sourceforge.net/project/swig/swig/swig-2.0.11/swig-2.0.11.tar.gz
tar -xvf swig-2.0.11.tar.gz
cd swig-2.0.11; ./configure --prefix=$LOCAL; make -j; make install; cd ..
ln -s 

pip install numpy >= 1.6.4
SWIG_FEATURES="-cpperraswarn -includeall -I/usr/include/openssl" pip install M2Crypto

# Install metaio
wget https://www.lsc-group.phys.uwm.edu/daswg/download/software/source/metaio-8.2.tar.gz
tar -xvf metaio-8.2.tar.gz
cd metaio-8.2; CPPFLAGS=-std=gnu99 ./configure --prefix=$LOCAL; make -j; make install; cd ..

# install framel
wget http://lappweb.in2p3.fr/virgo/FrameL/v8r26.tar.gz
tar -xvf v8r26.tar.gz 
cd v8r26; autoreconf; ./configure --prefix=$LOCAL;make -j; make install; cd ..

cd usr/lib/pkgconfig
sed -i "s|usr|$LOCAL/usr|g" gsl.pc

# Install lalsuite itself
cd $LOCAL/src/
git clone https://github.com/ahnitz/lalsuite.git
cd lalsuite
./00boot; ./configure --prefix=$LOCAL --enable-swig-python; make -j; make install
source $LOCAL/etc/lal-user-env.sh

echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> $TRAVIS_BUILD_DIR/source
echo "export PKG_CONFIG_PATH=$PKG_CONFIG_PATH" >> $TRAVIS_BUILD_DIR/source
echo "export PATH=$PATH" >> $TRAVIS_BUILD_DIR/source
echo "export HDF5_DIR=$HDF5_DIR" >> $TRAVIS_BUILD_DIR/source
echo "export LAPACK=$LOCAL/usr/lib/lapack" >> $TRAVIS_BUILD_DIR/source
echo "export BLAS=$LOCAL/usr/lib/libblas" >> $TRAVIS_BUILD_DIR/source
echo source $LOCAL/etc/lal-user-env.sh >> $TRAVIS_BUILD_DIR/source
echo $PWD
chmod 755 $TRAVIS_BUILD_DIR/source
