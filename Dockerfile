FROM python:3.10-slim-bullseye


RUN apt-get update && \
    apt-get install -y build-essential && \
    apt-get clean

## to build from local sources, use the lines below:
# RUN pip install numba fast-fisher pandas scipy scikit-learn matplotlib statsmodels fire
# WORKDIR /tmp/scoary
# COPY . .
# RUN pip install .
# RUN rm -rf /tmp/scoary
# WORKDIR /

ARG SCOARY_VERSION
RUN pip install scoary-2==$SCOARY_VERSION && \
    pip cache purge

# set these environment variables to directories where non-root is allowed to write
ENV NUMBA_CACHE_DIR=/tmp/NUMBA_CACHE_DIR
ENV CONFINT_DB=/tmp/CONFINT_DB
ENV MPLCONFIGDIR=/tmp/MPLCONFIGDIR
