FROM nfcore/base:1.7

LABEL description="Docker image containing all requirements for nf-core/eager pipeline"
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-eager-2.1.0dev/bin:$PATH
