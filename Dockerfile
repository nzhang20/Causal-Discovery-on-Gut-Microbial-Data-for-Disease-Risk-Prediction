# 1) choose base container
# generally use the most recent tag

# base notebook, contains Jupyter and relevant tools
# See https://github.com/ucsd-ets/datahub-docker-stack/wiki/Stable-Tag 
# for a list of the most current containers we maintain
ARG BASE_CONTAINER=ghcr.io/ucsd-ets/datascience-notebook:stable

FROM $BASE_CONTAINER

LABEL maintainer="UC San Diego ITS/ETS <ets-consult@ucsd.edu>"

# 2) change to root to install packages
USER root

RUN apt update

RUN apt-get -y install traceroute nmap aria2

# 3) install packages using notebook user
USER jovyan

RUN conda install -y python-graphviz conda-forge::cmdstan

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

RUN R -e "install.packages( \
	c('glasso', 'dplyr', 'data.table', 'rjson', 'glmnet'), \
	repos='http://cran.rstudio.com/')"

# Override command to disable running jupyter notebook at launch
CMD ["/bin/bash"]
