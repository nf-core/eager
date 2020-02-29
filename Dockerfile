FROM nfcore/base:1.9
LABEL authors="James A. Fellows Yates, Thiseas C. Lamnidis, Maxime Borry, Maxime Garcia, Aida Andrades Valtuena, Zandra Fagernäs, Stephen Clayton, Judith Neukamm, Alexander Peltzer" \
      description="Docker image containing all requirements for nf-core/eager pipeline"
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-eager-2.1.0dev > nf-core-eager-2.1.0dev.yml
ENV PATH /opt/conda/envs/nf-core-eager-2.1.0dev/bin:$PATH
