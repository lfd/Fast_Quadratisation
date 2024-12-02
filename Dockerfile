# Start your image with a node base image
FROM ubuntu:22.04

LABEL maintainer="repropackage@qsw24.github.com"

# The /repro directory should act as the main application directory
WORKDIR /repro



#Update
RUN apt-get update \
    && apt-get install -y build-essential \
    && apt-get install -y wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

#Install conda
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda
	
ENV PATH=$CONDA_DIR/bin:$PATH

RUN apt update && apt upgrade
RUN apt-get install -y apt-transport-https
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y libssl-dev
RUN apt-get install -y libfontconfig1-dev
RUN apt-get install -y libharfbuzz-dev libfribidi-dev
RUN apt-get install -y libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
#RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/r-project.gpg
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" >> /etc/apt/sources.list 
RUN apt-get update
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get install -y r-base

RUN R -e "install.packages('tidyverse', dependencies=TRUE)"
RUN R -e "install.packages('patchwork', dependencies=TRUE)"
RUN R -e "install.packages('tikzDevice', dependencies=TRUE)"
RUN R -e "install.packages('scales', dependencies=TRUE)"
RUN R -e "install.packages('ggh4x', dependencies=TRUE)"
RUN R -e "install.packages('ggpmisc', dependencies=TRUE)"

RUN apt-get install -y texlive-full

RUN conda create -n quark_install python=3.10

# Make RUN commands use quark_install:
SHELL ["conda", "run", "-n", "quark_install", "/bin/bash", "-c"]

# git clone https://gitlab.com/quantum-computing-software/quark.git
# cd quark
# conda create -n quark_install python=3.10
# conda activate quark_install
# conda config --add channels conda-forge
# conda install --file ./requirements.txt
# conda install conda-build
# conda develop ./

RUN conda config --add channels conda-forge
RUN conda update -n base -c defaults conda
RUN conda install -c dlr-sc quark=1.1
RUN conda install -c conda-forge pathos
RUN conda install pandas
RUN conda install scipy
RUN pip install p_tqdm
RUN pip install sortedcontainers

# Copy local directories to the current local directory of our docker image (/repro)
COPY ./code/ ./

#OpenShell
CMD ["/bin/bash"]

#CMD ["source", "activate", "quark_install"]
