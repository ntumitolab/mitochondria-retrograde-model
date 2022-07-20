FROM ghcr.io/sosiristseng/docker-juliabook:1.8.0rc3

WORKDIR /usr/src/app

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# Julia dependencies
COPY Project.toml Manifest.toml ./
COPY src/ src
RUN julia --threads=auto --color=yes --project=@. -e 'import Pkg; Pkg.instantiate()'

CMD ["jupyter-book"]
