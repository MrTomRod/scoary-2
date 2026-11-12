FROM ghcr.io/astral-sh/uv:python3.14-trixie-slim

ARG SCOARY_VERSION
ENV UV_COMPILE_BYTECODE=1 UV_NO_CACHE=1

# to build from local sources (run `uv build` first), use the lines below:
COPY dist/scoary_2-$SCOARY_VERSION-py3-none-any.whl /tmp/scoary/
RUN uv pip install --system /tmp/scoary/scoary_2-$SCOARY_VERSION-py3-none-any.whl && \
    rm -rf /tmp/scoary

# to build from PyPI, use this instead:
# RUN uv pip install --system scoary-2==$SCOARY_VERSION

# set these environment variables to directories where non-root is allowed to write
ENV NUMBA_CACHE_DIR=/tmp/NUMBA_CACHE_DIR
ENV CONFINT_DB=/tmp/CONFINT_DB
ENV MPLCONFIGDIR=/tmp/MPLCONFIGDIR

WORKDIR /data
