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
# 1. System dependencies (all libraries needed by any R package in this image)
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl wget ca-certificates \
        perl \
        ncbi-blast+ \
        python3 python3-pip python3-dev \
        build-essential pkg-config \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        libwebp-dev \
        libhdf5-dev \
        libgit2-dev \
        libssh2-1-dev \
        libglpk-dev \
        libsodium-dev \
        pandoc \
        git \
        emboss \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 2. IgBLAST + bundled IMGT germlines (from Immcantation suite image)
# ---------------------------------------------------------------------------
ARG IGBLAST_VERSION=1.22.0
RUN curl -fsSL \
    "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz" \
    | tar -xz -C /opt \
    && ln -s /opt/ncbi-igblast-${IGBLAST_VERSION} /opt/igblast \
    && ln -s /opt/igblast/bin/* /usr/local/bin/

ENV IGDATA=/opt/igblast

COPY --from=immcantation/suite:4.5.0 \
     /usr/local/share/germlines \
     /opt/igblast/germlines

ENV IGBLAST_GERMLINES=/opt/igblast/germlines

# ---------------------------------------------------------------------------
# 3. Immcantation Python tools (Change-O / MakeDb.py)
# ---------------------------------------------------------------------------
RUN pip3 install --no-cache-dir changeo presto

# ---------------------------------------------------------------------------
# 4. R: configure .libPaths() to include previous R version libraries.
#    This lets R find packages already installed by the base rocker image
#    (e.g. from R 4.3.x site-library) without reinstalling them, and makes
#    the build resilient to minor R version bumps.
#    Written to Rprofile.site so every Rscript call in this image sees it.
# ---------------------------------------------------------------------------
RUN Rscript -e "  r_ver  <- paste(R.version\$major, R.version\$minor, sep='.');   r_maj  <- R.version\$major;   # Candidate previous-version library paths (site and user variants)
  cands  <- c(     sprintf('/usr/local/lib/R/library'),                        sprintf('/usr/lib/R/library'),                              sprintf('/usr/local/lib/R/site-library'),                   sprintf('/home/%s/R/x86_64-pc-linux-gnu-library', Sys.getenv('USER')),     sprintf('/root/R/x86_64-pc-linux-gnu-library/%s', r_ver),     sprintf('/usr/local/lib/R/site-library')                  );   existing <- unique(c(.libPaths(), cands[dir.exists(cands)]));   lib_line <- paste0('.libPaths(c(',     paste0('"', existing, '"', collapse=','),     '))');   rprofile <- file.path(R.home('etc'), 'Rprofile.site');   write(lib_line, rprofile, append=TRUE);   message('Rprofile.site updated with ', length(existing), ' library paths:');   message(paste(' ', existing, collapse='
'))"

# ---------------------------------------------------------------------------
# 5. R: install pak — resolves the full dependency graph automatically,
#    installs packages in the correct order, and handles binary vs source.
# ---------------------------------------------------------------------------
RUN Rscript -e "  install.packages('pak', repos='https://cloud.r-project.org');   stopifnot(requireNamespace('pak', quietly=TRUE));   message('pak ', packageVersion('pak'), ' ready')"

# ---------------------------------------------------------------------------
# 6. R: Bioconductor dependencies via pak
# ---------------------------------------------------------------------------
RUN Rscript -e "\
  pak::pkg_install(c( \
    'bioc::Biostrings','bioc::DECIPHER', \
    'bioc::GenomicAlignments','bioc::GenomicRanges', \
    'bioc::IRanges','bioc::S4Vectors', \
    'bioc::ComplexHeatmap' \
  ), ask=FALSE)"

# ---------------------------------------------------------------------------
# 7. R: devtools and alakazam/tigger via pak
#    pak resolves the entire graph (fs, bslib, ragg, etc.) automatically
# ---------------------------------------------------------------------------
RUN Rscript -e "\
  pak::pkg_install(c('devtools','alakazam','tigger'), ask=FALSE); \
  stopifnot(requireNamespace('devtools', quietly=TRUE)); \
  message('devtools ', packageVersion('devtools'), ' ready')"

# ---------------------------------------------------------------------------
# 8. R: remaining blendAIRR runtime dependencies
# ---------------------------------------------------------------------------
RUN Rscript -e "\
  pak::pkg_install(c( \
    'plotly','dplyr','ggplot2','R6','jsonlite','magrittr','dendextend', \
    'RColorBrewer','circlize','splitstackshape','zen4R','rlang', \
    'data.table','optparse','remotes', \
    'Rcpp','stringdist','igraph','cluster','ape' \
  ), ask=FALSE)"

# ---------------------------------------------------------------------------
# 9. PIgLET — copied directly from a pre-built local R library.
#
#    This avoids all compilation and Bitbucket auth issues.
#    The compiled package (including piglet.so) is copied verbatim into
#    the image's R site-library.
#
#    How to prepare vendor/piglet-built/:
#      # On your local machine (where piglet is already installed):
#      Rscript -e ".libPaths()"           # note your R library path
#      cp -r <R_LIB_PATH>/piglet vendor/piglet-built
#      git add vendor/piglet-built
#      git commit -m "vendor pre-built piglet from local R library"
#
#    To update piglet: reinstall locally, repeat the cp above.
# ---------------------------------------------------------------------------
COPY vendor/piglet-built /usr/local/lib/R/site-library/piglet

RUN Rscript -e "\
  library(piglet); \
  message('PIgLET OK: ', packageVersion('piglet'))"

# ---------------------------------------------------------------------------
# 10. Verify all required packages load
# ---------------------------------------------------------------------------
RUN Rscript -e "\
  pkgs <- c('piglet','Biostrings','DECIPHER','data.table','optparse', \
            'alakazam','tigger','devtools'); \
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; \
  if (length(missing)) stop('Missing: ', paste(missing, collapse=', ')); \
  message('All R packages OK')"

# ---------------------------------------------------------------------------
# 11. Install blendAIRR scripts
# ---------------------------------------------------------------------------
COPY build_hybrid_igblast_ref.sh   /usr/local/bin/build_hybrid_igblast_ref
COPY R/piglet_annotate_and_build.R /opt/blendAIRR/R/piglet_annotate_and_build.R

RUN chmod +x /usr/local/bin/build_hybrid_igblast_ref

ENV HYBRID_IGBLAST_R_SCRIPT=/opt/blendAIRR/R/piglet_annotate_and_build.R

# ---------------------------------------------------------------------------
# 12. Final smoke-test
# ---------------------------------------------------------------------------
RUN igblastn -version && \
    Rscript -e "library(piglet); cat('piglet', as.character(packageVersion('piglet')), '\n')" && \
    ls /opt/igblast/germlines/imgt/ | head -3

# ---------------------------------------------------------------------------
# 13. Entry point
# ---------------------------------------------------------------------------
WORKDIR /data

ENTRYPOINT ["build_hybrid_igblast_ref"]
CMD ["--help"]
