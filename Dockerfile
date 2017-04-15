FROM andrewosh/binder-base

USER root
RUN echo 'deb http://software.ligo.org/lscsoft/debian/ jessie contrib' >> /etc/apt/sources.list
RUN apt-get update
RUN apt-get --force-yes install lalsimulation

USER main

RUN pip install pycbc --user
