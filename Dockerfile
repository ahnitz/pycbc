FROM andrewosh/binder-base

USER root

RUN apt-get update
RUN echo 'deb http://software.ligo.org/lscsoft/debian/ jessy universe' >> /etc/apt/sources.list
RUN apt-get install lscsoft-all
