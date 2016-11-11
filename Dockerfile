#
#  From this base-image / starting-point
#
FROM debian:testing

#
#  Authorship
#
MAINTAINER ap13@sanger.ac.uk

#
# A shared directory for results
#
VOLUME ["/data"]

#
# Update and Install dependencies
#
RUN apt-get update -qq && apt-get install -y python-genometools genometools wget python-dev python-setuptools

#
# Download build and install gff3toembl python
#
RUN wget https://github.com/sanger-pathogens/gff3toembl/archive/v1.1.4.tar.gz && mv v1.1.4.tar.gz /opt && cd /opt && tar xzf v1.1.4.tar.gz

RUN cd /opt && tar xzf v1.1.4.tar.gz && rm v1.1.4.tar.gz && cd /opt/gff3toembl-1.1.4 && python setup.py install
