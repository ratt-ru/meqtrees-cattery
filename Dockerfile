FROM kernsuite/base:2
RUN docker-apt-install \
    casacore-dev \
    casacore-tools \
    casarest \
    python-pip \
    python-meqtrees-timba \
    python-pyfits \
    python-numpy \
    python-scipy \
    lofar \
    python-astlib \
    python-purr \
    meqtrees-timba \
    python-kittens \
    python-owlcat \
    python-tigger
ADD . /code
WORKDIR /code
RUN python setup.py install
ENV MEQTREES_CATTERY_PATH=/usr/local/lib/python2.7/dist-packages/Cattery
WORKDIR /code/test/Batchtest
RUN python batch_test.py
