FROM nfcore/base:1.9
LABEL authors="Alexander Peltzer, Stephen Clayton, James A. Fellows Yates, Maxime Borry, Thiseas C. Lamnidis, Zandra Fagernäs" \
      description="Docker image containing all software requirements for the nf-core/eager pipeline"

LABEL description="Docker image containing all requirements for nf-core/eager pipeline"
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-eager-2.1.0dev > nf-core-eager-2.1.0dev.yml
ENV PATH /opt/conda/envs/nf-core-eager-2.1.0dev/bin:$PATH
