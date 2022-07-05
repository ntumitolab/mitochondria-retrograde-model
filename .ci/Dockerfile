FROM python:3.10-slim

WORKDIR /usr/src/app

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# System
RUN apt-get update && apt-get install -y build-essential git gzip brotli parallel --no-install-recommends && rm -rf /var/lib/apt/lists/*

# Julia
ENV JULIA_PATH /usr/local/julia/
ENV PATH $JULIA_PATH/bin:$PATH
COPY --from=julia:1.8.0-rc1 ${JULIA_PATH} /usr/local/

# Python dependencies
COPY .ci/requirements.txt .ci/requirements.txt
RUN pip install --no-cache-dir -U pip wheel setuptools \
    && pip install --no-cache-dir -r .ci/requirements.txt

# Julia dependencies
COPY Project.toml Manifest.toml ./
COPY src/ src
RUN julia --threads=auto --color=yes --project="" -e 'import Pkg; Pkg.add("IJulia")' && julia --threads=auto --color=yes --project=@. -e 'import Pkg; Pkg.instantiate()'

CMD ["jupyter-book"]
