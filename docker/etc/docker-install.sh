#!/bin/bash -v
set -e
cd /scratch
python3.10 -m pip install --upgrade pip
python3.10 -m pip install -r requirements.txt
python3.10 -m pip install -r requirements-igwn.txt
python3.10 -m pip install -r companion.txt
python3.10 -m pip install .
cd /
mkdir -p /opt/pycbc/src
cp -a /scratch /opt/pycbc/src/pycbc
chmod a+rw /dev/fuse
mount /cvmfs/config-osg.opensciencegrid.org
mount /cvmfs/software.igwn.org
mkdir -p /opt/pycbc/pycbc-software/share/lal-data
rsync --exclude='SEOBNRv1ROM*' --exclude='SEOBNRv2ROM_DS_HI_v1.0.hdf5' --exclude='NRSur4d2s_FDROM.hdf5' -ravz /cvmfs/software.igwn.org/pycbc/lalsuite-extra/current/share/lalsimulation/ /opt/pycbc/pycbc-software/share/lal-data/
umount /cvmfs/config-osg.opensciencegrid.org
umount /cvmfs/software.igwn.org
chown -R 1000:1000 /opt/pycbc/src /opt/pycbc/pycbc-software
chmod -R u=rwX,g=rX,o=rX /opt/pycbc
exit 0
