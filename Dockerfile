FROM python:3.10-slim-bullseye


RUN apt-get update && \
    apt-get install -y build-essential && \
    apt-get clean

ARG SCOARY_VERSION

# to build from local sources, use the lines below:
COPY dist/*$SCOARY_VERSION* /tmp/scoary/
RUN pip install -U /tmp/scoary/scoary_2-$SCOARY_VERSION-py3-none-any.whl && \
    pip cache purge && \
    rm -rf /tmp/scoary

# to build from pip, use this:
# RUN pip install scoary-2==$SCOARY_VERSION && \
#     pip cache purge

# set these environment variables to directories where non-root is allowed to write
ENV NUMBA_CACHE_DIR=/tmp/NUMBA_CACHE_DIR
ENV CONFINT_DB=/tmp/CONFINT_DB
ENV MPLCONFIGDIR=/tmp/MPLCONFIGDIR
