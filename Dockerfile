FROM ubuntu:22.04

LABEL maintainer="lukas.schmidbauer@othr.com"

WORKDIR /repro

RUN apt-get -y update \
    && apt-get install -y build-essential \
    && apt-get install -y wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get -y update && apt-get -y upgrade
RUN apt-get install -y apt-transport-https libcurl4-openssl-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libwebp-dev

ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2025.12-2-Linux-x86_64.sh -O ~/conda.sh  && \
    bash ~/conda.sh -b -p /opt/conda
	
ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" >> /etc/apt/sources.list 
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get install -y r-base cmake

RUN R -e "install.packages('tidyverse', dependencies=TRUE)"
RUN R -e "install.packages('patchwork', dependencies=TRUE)"
RUN R -e "install.packages('tikzDevice', dependencies=TRUE)"
RUN R -e "install.packages('scales', dependencies=TRUE)"
RUN R -e "install.packages('ggh4x', dependencies=TRUE)"

RUN apt-get install -y texlive-full p7zip-full

RUN conda create -n quark_install python=3.10

SHELL ["conda", "run", "-n", "quark_install", "/bin/bash", "-c"]

RUN conda config --add channels conda-forge
RUN conda update -n base -c defaults conda -vvv
RUN conda install -c dlr-sc quark=1.1
RUN conda install -c conda-forge pathos
RUN conda install conda-forge::scip==9.0.0
RUN pip install pandas==2.2.2
RUN pip install scipy==1.13.1
RUN pip install p_tqdm==1.4.0
RUN pip install sortedcontainers==2.4.0
RUN pip install qiskit[visualization]==2.1.0
RUN pip install multiset==3.2.0
RUN pip install networkx==3.3

COPY ./code/ ./

CMD ["/bin/bash"]