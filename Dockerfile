FROM continuumio/miniconda3

RUN conda install gcc_linux-64 gxx_linux-64 python=3.10
RUN ln -s /opt/conda/bin/x86_64-conda_cos7-linux-gnu-gcc /usr/bin/gcc

ARG SCOARY_VERSION
RUN pip install scoary-2==$SCOARY_VERSION && \
    pip cache purge

# set these environment variables to directories where non-root is allowed to write
ENV NUMBA_CACHE_DIR=/tmp/NUMBA_CACHE_DIR
ENV CONFINT_DB=/tmp/CONFINT_DB
ENV MPLCONFIGDIR=/tmp/MPLCONFIGDIR

WORKDIR /data
