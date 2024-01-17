FROM python:3.11-bookworm

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

RUN pip install git+https://github.com/cerebis/sim3C.git

RUN mkdir -p /app
WORKDIR /app
ENTRYPOINT ["sim3C"]
CMD ["--help"]
