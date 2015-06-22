#!/bin/bash
#set -e

LOCAL=$PWD
mkdir -p $LOCAL/src

export LD_LIBRARY_PATH=$LOCAL/lib:$LOCAL/lib64:$LOCAL/usr/lib:$LOCAL/usr/lib64
export PKG_CONFIG_PATH=$LOCAL/lib/pkgconfig
export PATH=$PATH:$LOCAL/bin:$LOCAL/usr/bin

apt-get download \
libfftw3-3 \
libfftw3-dev \
libhdf5-serial-1.8.4 \
libhdf5-serial-dev \
libblas3gf \
libblas-dev \
liblapack3gf \
liblapack-dev \
gfortran \
libgsl0ldbl \
libgsl0-dev

for PKG in *.deb; do
    dpkg-deb -x $PKG $LOCAL
done

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

# Install lalsuite itself
cd $LOCAL/src/
git clone https://github.com/ahnitz/lalsuite.git
cd lalsuite
./00boot; ./configure --prefix=$LOCAL --enable-swig-python; make -j; make install
source $LOCAL/etc/lal-user-env.sh

echo source $LOCAL/etc/glue-user-env.sh >> $TRAVIS_BUILD_DIR/source
echo source $LOCAL/etc/lal-user-env.sh >> $TRAVIS_BUILD_DIR/source
echo $PWD
chmod 755 $TRAVIS_BUILD_DIR/source
