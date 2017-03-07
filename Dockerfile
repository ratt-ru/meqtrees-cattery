FROM kernsuite/base:dev
RUN docker-apt-install \
    casacore-dev \
    casacore-tools \
    casarest \
    python-pip \
    python-meqtrees-timba \
    python-pyfits \
    python-numpy \
    python-scipy \
    python-astlib \
    python-purr \
    meqtrees-timba \
    python-kittens \
    python-owlcat \
    python-tigger \
    python-meqtrees-cattery \
    python-pyxis
RUN docker-apt-install \
    makems
ADD . /code
WORKDIR /code
RUN python setup.py install
WORKDIR /code/test/Batchtest
CMD python batch_test.py
