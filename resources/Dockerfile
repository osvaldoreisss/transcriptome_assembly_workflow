FROM ubuntu:latest
LABEL Author: Osvaldo Reis Junior <osvaldoreisss@unicamp.br>
WORKDIR /opt
ENV PATH /opt/conda/bin:${PATH}

RUN apt-get update && apt-get install -y build-essential
RUN apt-get install -y git wget

ADD https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh .

RUN bash Miniconda2-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda2-latest-Linux-x86_64.sh

RUN git clone https://github.com/harta55/EnTAP.git

# Install Genemark
ADD http://topaz.gatech.edu/GeneMark/tmp/GMtool_CLFcF/gm_key_64.gz .
ADD http://topaz.gatech.edu/GeneMark/tmp/GMtool_CLFcF/gmst_linux_64.tar.gz .
RUN gunzip gm_key_64.gz && mv gm_key_64 /root
RUN mkdir EnTAP/libs/gmst_linux_64/ && tar -zxvf  gmst_linux_64.tar.gz --directory EnTAP/libs/gmst_linux_64/
ENV PATH /opt/EnTAP/libs/gmst_linux_64/:${PATH}

# Install Interproscan
RUN mkdir my_interproscan && cd my_interproscan && \
    wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.47-82.0/interproscan-5.47-82.0-64-bit.tar.gz && \
    wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.47-82.0/interproscan-5.47-82.0-64-bit.tar.gz.md5 && \
    md5sum -c interproscan-5.47-82.0-64-bit.tar.gz.md5 && \
    tar -pxvzf interproscan-5.47-82.0-*-bit.tar.gz && \
    cd interproscan-5.47-82.0 && python3 initial_setup.py

# Install Diamond
RUN cd EnTAP/libs/diamond-0.8.31 && \
    mkdir bin && \
    cd bin && cmake .. && make && make install

# Install RSEM
RUN cd EnTAP/libs/RSEM-1.3.0 && \
    make && make ebseq && make install

# Install Entap
RUN cd EnTAP && cmake CMakeLists.txt -DCMAKE_INSTALL_PREFIX=/usr/bin && \
    make && make install
