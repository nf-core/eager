From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Alexander Peltzer <alexander.peltzer@qbic.uni-tuebingen.de>
    DESCRIPTION Container image containing all requirements for the nf-core/EAGER2 pipeline
    VERSION 2.0dev

%files
    environment.yml /

%post
    /opt/conda/bin/conda env update -n root -f /environment.yml 
    /opt/conda/bin/conda clean -a
