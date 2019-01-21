FROM kernsuite/base:5
RUN docker-apt-install \
    casacore-dev \
    casacore-tools \
    casarest \
    python-pip \
    python-pyfits \
    python-numpy \
    python-scipy \
    python-astlib \
    python-casacore \
    libqdbm-dev \
    build-essential \
    cmake \
    libblitz0-dev \
    binutils-dev \
    makems

RUN docker-apt-install \
    wget

################################
# install latest masters
################################
RUN echo "deb-src http://ppa.launchpad.net/kernsuite/kern-5/ubuntu bionic main" > /etc/apt/sources.list.d/kernsuite-ubuntu-kern-4-bionic.list
RUN apt-get update
RUN apt-get build-dep -y meqtrees-timba meqtrees-cattery
RUN docker-apt-install python-qt4 python-qwt5-qt4 git
RUN mkdir -p /opt/src/meqtrees
ENV BUILD /opt/src/meqtrees
WORKDIR $BUILD

# now the rest of Meqtrees
WORKDIR $BUILD
RUN git clone https://github.com/ska-sa/meqtrees-timba
RUN git clone https://github.com/ska-sa/kittens 
RUN git clone https://github.com/ska-sa/purr
RUN git clone https://github.com/ska-sa/tigger 
RUN git clone https://github.com/ska-sa/owlcat 

# Add current module
ADD . /opt/src/meqtrees/cattery 

# Build Timba
RUN mkdir meqtrees-timba/build
WORKDIR $BUILD/meqtrees-timba
RUN Tools/Build/bootstrap_cmake release
WORKDIR $BUILD
RUN cd meqtrees-timba/build/release && make -j16 && make install

# Install python2.7 modules
RUN cd cattery && pip install . 
RUN cd kittens && pip install .
RUN cd purr && pip install .
RUN cd tigger && pip install .
RUN cd owlcat && pip install .

ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH
ENV PATH /usr/local/bin:$PATH
RUN ldconfig # update the linker config

################################
# get the test from pyxis
################################
WORKDIR $BUILD
RUN wget https://github.com/ska-sa/pyxis/archive/v1.6.2.tar.gz
RUN tar -xvf v1.6.2.tar.gz
WORKDIR $BUILD/pyxis-1.6.2
RUN ls .
RUN pip install -U .

# run test when built
RUN pip install nose
WORKDIR /usr/local/lib/python2.7/dist-packages/Pyxis/recipies/meqtrees-batch-test
RUN python2.7 -m "nose"

ENTRYPOINT ["meqtree-pipeliner.py"]
CMD ["--help"]