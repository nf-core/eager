From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Alexander Peltzer <alexander.peltzer@qbic.uni-tuebingen.de>
    DESCRIPTION Container image containing all requirements for the nf-core/eager pipeline
    VERSION 2.0.6

%environment
    PATH=/opt/conda/envs/nf-core-eager-2.0.6/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml 
    /opt/conda/bin/conda clean -a
