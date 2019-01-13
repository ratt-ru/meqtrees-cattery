FROM kernsuite/base:4
RUN docker-apt-install \
    casacore-dev \
    casacore-tools \
    casarest \
    python-pip \
    python3 \
    python3-pip \
    python-pyfits \
    python-numpy \
    python-scipy \
    python-astlib \
    python-casacore \
    libqdbm-dev
RUN docker-apt-install \
    build-essential \
    cmake \
    lofar

################################
# install latest masters
################################
RUN apt-get update
RUN apt-get build-dep -y meqtrees-timba meqtrees-cattery
RUN docker-apt-install python-qt4 python-qwt5-qt4 git
RUN mkdir -p /opt/src/meqtrees
ENV BUILD /opt/src/meqtrees
WORKDIR $BUILD
RUN git clone https://github.com/ska-sa/meqtrees-timba
RUN git clone https://github.com/ska-sa/kittens 
RUN git clone https://github.com/ska-sa/purr 
RUN git clone https://github.com/ska-sa/pyxis 
RUN git clone https://github.com/ska-sa/tigger 
RUN git clone https://github.com/ska-sa/meqtrees-cattery 
RUN git clone https://github.com/ska-sa/owlcat 
RUN mkdir meqtrees-timba/build
WORKDIR $BUILD/meqtrees-timba
RUN Tools/Build/bootstrap_cmake release
WORKDIR $BUILD
RUN cd meqtrees-timba/build/release && make -j8 && make install

#Install python3
RUN cd kittens && python3 setup.py install
RUN cd purr && python3 setup.py install
RUN cd pyxis && python3 setup.py install
RUN cd tigger && python3 setup.py install
RUN cd meqtrees-cattery && python3 setup.py install
RUN cd owlcat && python3 setup.py install

#Install python2.7
RUN cd kittens && python2.7 setup.py install
RUN cd purr && python2.7 setup.py install
RUN cd pyxis && python2.7 setup.py install
RUN cd tigger && python2.7 setup.py install
RUN cd meqtrees-cattery && python2.7 setup.py install
RUN cd owlcat && python2.7 setup.py install

ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH
ENV PATH /usr/local/bin:$PATH
RUN ldconfig # update the linker config
################################
# add test
################################
ADD . /code
WORKDIR /code/test/Batchtest
CMD echo '##################################################'
CMD echo 'Starting Python 2.7 testing...'
CMD echo '##################################################'
CMD python2.7 batch_test.py
CMD echo '##################################################'
CMD echo 'Starting Python 3 testing...'
CMD echo '##################################################'
CMD python3 batch_test.py
