FROM broadinstitute/gatk:4.2.0.0

LABEL maintainer="https://github.com/p4rkerw"

# update gpg key to enable apt-get update
RUN wget https://packages.cloud.google.com/apt/doc/apt-key.gpg && apt-key add apt-key.gpg

# install htslib and bcftools
ARG SAMTOOLS_VERSION=1.9
ARG BCFTOOLS_VERSION=1.9
ARG SHAPEIT_VERSION=v4.2.0

# update ubuntu packages
RUN apt update && apt install -y \
                 ca-certificates \
                 libbz2-dev \
                 libcurl4-openssl-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 autoconf \
                 automake \
                 bzip2 \
                 gcc \
		             g++ \
                 make \
                 wget \
                 zlib1g-dev \
		             libboost-all-dev \
		             git \
		             pv

# install htslib via samtools
RUN \                
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && ./configure CFLAGS="-static" --without-curses && make -s all all-htslib && make install install-htslib && cd - && \ 
    rm -rf samtools-${SAMTOOLS_VERSION}*

# clone shapeit4 repo and compile with htslib. There is a conda shapeit4 option, but it doesn't install correctly
RUN \
    cd /opt && \   
    git clone --depth 1 --branch ${SHAPEIT_VERSION} https://github.com/odelaneau/shapeit4.git && \
    cd shapeit4 && \
    cat makefile | sed 's|HTSLIB_INC=$(HOME)/Tools/htslib-1.9|HTSLIB_INC=/usr/local/include|' | sed 's|HTSLIB_LIB=$(HOME)/Tools/htslib-1.9/libhts.a|HTSLIB_LIB=/usr/local/lib/libhts.a|' > makefile2 && \
    mv makefile makefile.old && \
    mv makefile2 makefile && \
    make && \
    ln -s /opt/shapeit4/bin/shapeit4.2 /usr/bin/shapeit4.2

# clone the 10X genomics subset-bam repo into /opt
RUN \
  cd /opt && \
  wget https://github.com/10XGenomics/subset-bam/releases/download/v1.1.0/subset-bam_linux && \
  chmod +x subset-bam_linux && \    
  ln -s /opt/subset-bam_linux /usr/bin/subset-bam

# install r bioconda pysam bwa STAR hdf5 for WASP dependencies
# pysam hdf5 required for WASP snp2h5 install and compile
RUN conda config --add channels r && \
    conda config --add channels bioconda && \
    conda install bwa=0.7.17 && \
    conda install bcftools=1.9 && \
    conda config --add channels conda-forge && \
    conda install hdf5 && \
    conda install pysam=0.15.3

# clone the wasp v0.3.4 repo and compile snp2h5 in /opt
RUN \
  cd /opt && \
  git clone -b v0.3.4 https://github.com/bmvdgeijn/WASP.git && \
  cat WASP/snp2h5/Makefile | sed 's/$(HDF_INSTALL)/\/opt\/miniconda/g' > WASP/snp2h5/Makefile2 && \
  mv WASP/snp2h5/Makefile2 WASP/snp2h5/Makefile && \
  export PATH=/opt/miniconda/bin/lib:$PATH && \
  cd WASP/snp2h5 && \
  make && \
  ln -s /opt/WASP/snp2h5/snp2h5 /usr/bin/snp2h5     

# update the snp2h5 path for hdf5 libs
ENV LD_LIBRARY_PATH "/opt/miniconda/lib:$LD_LIBRARY_PATH"

# install separately from main conda layer o/w tables module wont install properly with pysam
# this may be a result of conda not finding hdf5
# this module is required to run WASP
RUN pip install tables

# install STAR compiled linux executable
RUN \
  cd /opt && \
  wget https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz && \
  tar -xzf 2.5.1b.tar.gz && \
  ln -s /opt/STAR-2.5.1b/bin/Linux_x86_64/STAR /usr/bin/STAR

# update entrypoint
WORKDIR $HOME
