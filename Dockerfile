FROM centos:centos7

RUN rpm -ivh http://software.ligo.org/lscsoft/scientific/7.2/x86_64/production/lscsoft-production-config-1.3-1.el7.noarch.rpm
RUN yum -y install lscsoft-backports-config
RUN yum -y install lscsoft-epel-config
RUN yum -y install lscsoft-ius-config
RUN yum -y install git2u-all lscsoft-all
RUN yum install -y zlib-devel libpng-devel libjpeg-devel libsqlite3-dev sqlite-devel db4-devel
RUN yum -y install tkinter libpng-devel lynx telnet
RUN yum -y install compat-glibc compat-glibc-headers
RUN yum -y install gd-devel audit-libs-devel libcap-devel nss-devel
RUN yum -y install xmlto asciidoc hmaccalc newt-devel 'perl(ExtUtils::Embed)' pesign elfutils-devel binutils-devel numactl-devel pciutils-devel
RUN yum -y install dejagnu sharutils gcc-gnat libgnat dblatex gmp-devel mpfr-devel libmpc-devel
RUN yum -y install libuuid-devel netpbm-progs nasm
RUN yum -y install libstdc++-static
RUN yum -y install gettext-devel avahi-devel dyninst-devel crash-devel latex2html emacs libvirt-devel
RUN yum -y install xmlto-tex patch
RUN yum -y install ant asciidoc xsltproc fop docbook-style-xsl.noarch
RUN yum -y install vim-enhanced

