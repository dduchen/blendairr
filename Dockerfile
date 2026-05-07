# =============================================================================
# blendAIRR — Docker image
# https://github.com/dduchen/blendairr
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
# 2. IgBLAST — pre-compiled binary + bundled germline databases
#    The tarball from NCBI includes:
#      bin/            igblastn, makeblastdb, edit_imgt_file.pl, ...
#      internal_data/  per-organism internal BLAST dbs (includes mouse)
#      optional_file/  J-gene aux files (mouse_gl.aux etc.)
#      database/       (empty; user provides their own search databases)
# ---------------------------------------------------------------------------
ARG IGBLAST_VERSION=1.22.0
RUN curl -fsSL \
    "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz" \
    | tar -xz -C /opt \
    && ln -s /opt/ncbi-igblast-${IGBLAST_VERSION} /opt/igblast \
    && ln -s /opt/igblast/bin/* /usr/local/bin/

# IGDATA must point to the igblast share that contains internal_data/ and optional_file/.
# The tarball already ships these for all built-in organisms (mouse, human, etc.).
ENV IGDATA=/opt/igblast

# Download IMGT reference germlines for all built-in organisms.
# These are the *reference* sequences blendAIRR co-clusters against —
# NOT the user's custom sequences. Stored inside the igblast share so
# the build script finds them at the standard path.
RUN mkdir -p /opt/igblast/germlines && \
    curl -fsSL \
    "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/database.tar.gz" \
    | tar -xz -C /opt/igblast/ || true
# Also fetch the IMGT germlines package shipped with igblast releases
RUN curl -fsSL \
    "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/germline_database.tar.gz" \
    | tar -xz -C /opt/igblast/ 2>/dev/null || \
    # Fallback: clone from the igblast companion repository
    (git clone --depth 1 \
      https://bitbucket.org/kleinstein/immcantation /tmp/immcantation 2>/dev/null && \
     cp -r /tmp/immcantation/igblast/germlines /opt/igblast/ 2>/dev/null && \
     rm -rf /tmp/immcantation) || true

# Download IMGT germlines directly (the canonical source blendAIRR expects).
# Layout: germlines/imgt/<organism>/vdj/  and  germlines/imgt/<organism>/leader/
RUN mkdir -p /opt/igblast/germlines/imgt && \
    curl -fsSL \
      "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/imgt_germlines.tar.gz" \
      | tar -xz -C /opt/igblast/germlines/ 2>/dev/null || true

# Immcantation ships a reliable copy of IMGT germlines in their Docker image.
# We copy just the germline directory from there as a definitive fallback.
COPY --from=immcantation/suite:4.5.0 \
     /usr/local/share/germlines \
     /opt/igblast/germlines

# ---------------------------------------------------------------------------
# 3. Immcantation (Change-O / MakeDb.py)
# ---------------------------------------------------------------------------
RUN pip3 install --no-cache-dir changeo presto

# ---------------------------------------------------------------------------
# 4. R packages — Bioconductor + CRAN
#    Split into separate RUN layers so a single package failure doesn't
#    invalidate the entire (expensive) R install cache.
# ---------------------------------------------------------------------------
RUN Rscript -e "\
  options(repos = c(CRAN = 'https://cloud.r-project.org')); \
  install.packages(c('BiocManager','remotes','data.table','optparse','tools'), \
                   repos='https://cloud.r-project.org', quiet=TRUE)"

RUN Rscript -e "\
  BiocManager::install(version='3.19', ask=FALSE, update=FALSE); \
  BiocManager::install(c('Biostrings','DECIPHER'), ask=FALSE, update=FALSE)"

# ---------------------------------------------------------------------------
# 5. PIgLET — from Bitbucket
#    Pin PIGLET_COMMIT to a specific commit hash for reproducibility.
#    Find hashes at: https://bitbucket.org/yaarilab/piglet/commits/
#    Leave as HEAD to always get the latest (less reproducible).
# ---------------------------------------------------------------------------
ARG PIGLET_COMMIT=HEAD
RUN Rscript -e "\
  message('Installing PIgLET from Bitbucket (ref=${PIGLET_COMMIT})...'); \
  remotes::install_bitbucket('yaarilab/piglet', ref='${PIGLET_COMMIT}', \
                              upgrade='never', quiet=FALSE); \
  library(piglet); \
  message('PIgLET installed and loads OK: ', packageVersion('piglet'))"

# ---------------------------------------------------------------------------
# 6. Verify all dependencies load before shipping the image
# ---------------------------------------------------------------------------
RUN Rscript -e "\
  pkgs <- c('piglet','Biostrings','DECIPHER','data.table','optparse'); \
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; \
  if (length(missing)) stop('Missing R packages: ', paste(missing, collapse=', ')); \
  message('All R packages OK')"

# ---------------------------------------------------------------------------
# 7. Install blendAIRR tool scripts
# ---------------------------------------------------------------------------
COPY build_hybrid_igblast_ref.sh   /usr/local/bin/build_hybrid_igblast_ref
COPY R/piglet_annotate_and_build.R /opt/blendAIRR/R/piglet_annotate_and_build.R

RUN chmod +x /usr/local/bin/build_hybrid_igblast_ref

ENV HYBRID_IGBLAST_R_SCRIPT=/opt/blendAIRR/R/piglet_annotate_and_build.R

# Point the build script to the reference germlines inside the image.
# blendAIRR will look for: $IGBLAST_GERMLINES/imgt/<species>/vdj/
ENV IGBLAST_GERMLINES=/opt/igblast/germlines

# ---------------------------------------------------------------------------
# 8. Final smoke-test: confirm igblastn + R + germlines are all accessible
# ---------------------------------------------------------------------------
RUN igblastn -version && \
    Rscript -e "library(piglet); cat('piglet', as.character(packageVersion('piglet')), '\n')" && \
    ls /opt/igblast/germlines/imgt/ | head -5

# ---------------------------------------------------------------------------
# 9. Entry point
# ---------------------------------------------------------------------------
WORKDIR /data

ENTRYPOINT ["build_hybrid_igblast_ref"]
CMD ["--help"]
