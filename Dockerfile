# =============================================================================
# blendAIRR — Docker image
#
# Builds a hybrid IgBLAST germline reference database by merging custom
# species germline sequences with the closest IMGT reference species.
# Novel alleles are jointly clustered with reference sequences using PIgLET
# to inherit correct IMGT gene-family designations.
#
# Layers:
#   1. rocker/r-ver    — pinned R version on Ubuntu
#   2. System deps     — BLAST+, IgBLAST, Perl, Python + Immcantation
#   3. R packages      — Bioconductor + PIgLET (from Bitbucket)
#   4. blendAIRR tools — build script and R annotation script
# =============================================================================

FROM rocker/r-ver:4.4.1

LABEL org.opencontainers.image.title="blendAIRR" \
      org.opencontainers.image.description="Build hybrid IgBLAST germline reference databases for custom or non-reference species" \
      org.opencontainers.image.source="https://github.com/dduchen/blendairr" \
      org.opencontainers.image.documentation="https://github.com/dduchen/blendairr/blob/main/README.md" \
      org.opencontainers.image.licenses="MIT"

# ---------------------------------------------------------------------------
# 1. System dependencies
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl wget ca-certificates \
        perl \
        ncbi-blast+ \
        python3 python3-pip python3-dev \
        build-essential libcurl4-openssl-dev libssl-dev libxml2-dev \
        libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        git \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 2. IgBLAST — pre-compiled binary from NCBI FTP
# ---------------------------------------------------------------------------
ARG IGBLAST_VERSION=1.22.0
RUN curl -fsSL \
    "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz" \
    | tar -xz -C /opt \
    && ln -s /opt/ncbi-igblast-${IGBLAST_VERSION} /opt/igblast \
    && ln -s /opt/igblast/bin/* /usr/local/bin/
ENV IGDATA=/opt/igblast

# ---------------------------------------------------------------------------
# 3. Immcantation (Change-O / MakeDb.py)
# ---------------------------------------------------------------------------
RUN pip3 install --no-cache-dir changeo==1.3.0 presto==0.7.2

# ---------------------------------------------------------------------------
# 4. R packages — Bioconductor + CRAN
# ---------------------------------------------------------------------------
RUN Rscript -e "\
  options(repos = c(CRAN = 'https://cloud.r-project.org')); \
  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); \
  BiocManager::install(version='3.19', ask=FALSE, update=FALSE); \
  BiocManager::install(c('Biostrings','DECIPHER'), ask=FALSE, update=FALSE); \
  install.packages(c('data.table','optparse','remotes','tools'), \
                   repos='https://cloud.r-project.org')"

# ---------------------------------------------------------------------------
# 5. PIgLET — from Bitbucket (pin PIGLET_COMMIT to a specific commit hash
#    for fully reproducible builds; find hashes at:
#    https://bitbucket.org/yaarilab/piglet/commits/ )
# ---------------------------------------------------------------------------
ARG PIGLET_COMMIT=HEAD
RUN Rscript -e "\
  remotes::install_bitbucket('yaarilab/piglet', ref='${PIGLET_COMMIT}', \
                              upgrade='never', quiet=FALSE)"

# ---------------------------------------------------------------------------
# 6. Install blendAIRR tool scripts
# ---------------------------------------------------------------------------
COPY build_hybrid_igblast_ref.sh  /usr/local/bin/build_hybrid_igblast_ref
COPY R/piglet_annotate_and_build.R /opt/blendAIRR/R/piglet_annotate_and_build.R

RUN chmod +x /usr/local/bin/build_hybrid_igblast_ref

# Explicit env var so the build script always finds the R file regardless
# of working directory inside the container.
ENV HYBRID_IGBLAST_R_SCRIPT=/opt/blendAIRR/R/piglet_annotate_and_build.R

# ---------------------------------------------------------------------------
# 7. Default working directory and entry point
# ---------------------------------------------------------------------------
WORKDIR /data

ENTRYPOINT ["build_hybrid_igblast_ref"]
CMD ["--help"]
