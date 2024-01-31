FROM python:3.11-bookworm as build

LABEL maintainer="matt.demaere@gmail.com"
LABEL org.label-schema.schema-version="1.0"
LABEL org.label-schema.name="cerebis/qc3c"
LABEL org.label-schema.description="sim3C: read-pair simulation of 3C-based sequencing methodologies (HiC, Meta3C, DNase-HiC)"
LABEL org.label-schema.url="http://github.com/cerebis/sim3C/"
LABEL org.label-schema.vcs-url="http://github.com/cerebis/sim3C/"
#LABEL org.label-schema.vcs-ref="07cc4cd7e938950f9242df255d365d065b62d0a8"
LABEL org.label-schema.version="0.5"
LABEL org.label-schema.docker.cmd="docker run -v /path/to/data:/app cerebis/sim3c --help"

RUN apt-get update &&  \
    apt-get upgrade -y && \
    apt-get install -y llvm && \
    pip install -U pip && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir git+https://github.com/cerebis/sim3C.git && \
    sim3C --help

FROM python:3.11-slim as run

COPY --from=build /usr/local/lib/python3.11/site-packages /usr/local/lib/python3.11/site-packages
COPY --from=build /usr/local/bin/sim3C /usr/local/bin/sim3C

RUN mkdir -p /app
WORKDIR /app
ENTRYPOINT ["sim3C"]
CMD ["--help"]
