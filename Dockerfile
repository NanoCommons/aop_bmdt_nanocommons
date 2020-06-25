FROM rstudio/r-base:3.6.3-centos7 

USER root

RUN yum check-update
RUN yum -y install openssl-devel

RUN mkdir /home/scripts
ADD dependencies.R /home/scripts/dependencies.R

WORKDIR /home/scripts/
RUN Rscript /home/scripts/dependencies.R

CMD R
