ARG GET_VERSION=20190315_patched
ARG ROOT_VERSION=6.08
ARG UBUNTU_VERSION=xenial

FROM elitpc/get:${GET_VERSION}-${UBUNTU_VERSION}-${ROOT_VERSION}
ADD . /opt/soft/TPCReco-src
WORKDIR /opt/soft/TPCReco-src
RUN apt-get update -qq \
    && apt-get -y install $(cat requirements_apt.txt) \
    less vim emacs nano \
    && rm -rf /var/lib/apt/lists/*
RUN yes | pip3 install --no-cache-dir -r requirements_pip3.txt
WORKDIR /opt/soft/TPCReco-src/build
RUN cmake -DCMAKE_INSTALL_PREFIX=/opt/soft/TPCReco ..
RUN make -j && make install && make clean
ENV PATH=${PATH}:/opt/soft/TPCReco/bin
WORKDIR /opt/soft/TPCReco/resources
