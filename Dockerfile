FROM nfcore/base
MAINTAINER Alexander Peltzer <alex.peltzer@gmail.com>
LABEL authors="alex.peltzer@gmail.com" \
    description="Docker image containing all requirements for nf-core/EAGER2 pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-eager-2.0dev/bin:$PATH
