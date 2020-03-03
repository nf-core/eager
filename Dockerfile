FROM nfcore/base:1.9
LABEL authors="The nf-core/eager community" \
      description="Docker image containing all requirements for nf-core/eager pipeline"
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-eager-2.1.0 > nf-core-eager-2.1.0.yml
ENV PATH /opt/conda/envs/nf-core-eager-2.1.0/bin:$PATH
