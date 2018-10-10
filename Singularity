From:nfcore/base
Bootstrap:docker

%labels
<<<<<<< HEAD
    MAINTAINER Alexander Peltzer <alexander.peltzer@qbic.uni-tuebingen.de>
    DESCRIPTION Container image containing all requirements for the nf-core/EAGER2 pipeline
    VERSION 2.0dev

%environment
    PATH=/opt/conda/envs/nf-core-eager-2.0dev/bin:$PATH
=======
    DESCRIPTION Singularity image containing all requirements for the nf-core/eager pipeline
    VERSION 1.0dev

%environment
    PATH=/opt/conda/envs/nf-core-eager-1.0dev/bin:$PATH
>>>>>>> TEMPLATE
    export PATH

%files
    environment.yml /

%post
<<<<<<< HEAD
    /opt/conda/bin/conda env create -f /environment.yml 
=======
    /opt/conda/bin/conda env create -f /environment.yml
>>>>>>> TEMPLATE
    /opt/conda/bin/conda clean -a
