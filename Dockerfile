FROM andrewosh/binder-base

RUN echo 'deb http://software.ligo.org/lscsoft/debian/ jessy universe' >> /etc/apt/sources.list
RUN apt-get install lscsoft-all
