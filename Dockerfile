FROM andrewosh/binder-base

RUN add-apt-repository "deb http://software.ligo.org/lscsoft/debian/ universe"
RUN apt-get install lscsoft-all
