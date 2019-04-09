FROM nfcore/base
LABEL authors="Alexander Peltzer, Stephen Clayton, James A Fellows-Yates" \
      description="Docker image containing all requirements for nf-core/eager pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-eager-2.0.6/bin:$PATH
