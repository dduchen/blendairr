# Makefile — blendAIRR local development targets
# Usage: make <target>

IMAGE        ?= blendairr
TAG          ?= local
REGISTRY     ?= ghcr.io/YOUR_ORG
DATA_DIR     ?= $(PWD)/data
IGBLAST_VER  ?= 1.22.0
PIGLET_COMMIT ?= HEAD

.PHONY: build test run push singularity clean help

help:
	@echo ""
	@echo "  blendAIRR — make targets"
	@echo ""
	@echo "  build          Build the Docker image locally ($(IMAGE):$(TAG))"
	@echo "  test           Smoke-test the local image (R packages, igblastn, MakeDb.py)"
	@echo "  run            Run blendAIRR (set DATA_DIR and ARGS)"
	@echo "  push           Tag and push to $(REGISTRY)"
	@echo "  singularity    Convert local Docker image to a Singularity .sif file"
	@echo "  clean          Remove local image"
	@echo ""
	@echo "  Variables (override on command line):"
	@echo "    DATA_DIR=$(DATA_DIR)"
	@echo "    TAG=$(TAG)"
	@echo "    REGISTRY=$(REGISTRY)"
	@echo "    PIGLET_COMMIT=$(PIGLET_COMMIT)   (pin to a Bitbucket commit hash)"
	@echo ""
	@echo "  Example:"
	@echo "    make run ARGS='--species mouse --input_dir /data/MRL --outdir /data/out'"
	@echo ""

build:
	docker build \
	  --build-arg IGBLAST_VERSION=$(IGBLAST_VER) \
	  --build-arg PIGLET_COMMIT=$(PIGLET_COMMIT) \
	  -t $(IMAGE):$(TAG) .

test: build
	@echo "--- Testing --help ---"
	docker run --rm $(IMAGE):$(TAG) --help
	@echo "--- Testing R packages ---"
	docker run --rm --entrypoint Rscript $(IMAGE):$(TAG) \
	  -e "library(piglet); library(DECIPHER); library(data.table); cat('R packages OK\n')"
	@echo "--- Testing igblastn ---"
	docker run --rm --entrypoint bash $(IMAGE):$(TAG) \
	  -c "which igblastn && igblastn -version && which makeblastdb"
	@echo "--- Testing MakeDb.py ---"
	docker run --rm --entrypoint bash $(IMAGE):$(TAG) -c "MakeDb.py --version"
	@echo "All tests passed."

run: build
	docker run --rm \
	  -v "$(DATA_DIR)":/data \
	  $(IMAGE):$(TAG) $(ARGS)

push: build
	docker tag $(IMAGE):$(TAG) $(REGISTRY)/$(IMAGE):$(TAG)
	docker push $(REGISTRY)/$(IMAGE):$(TAG)

singularity: build
	singularity build $(IMAGE)_$(TAG).sif docker-daemon://$(IMAGE):$(TAG)
	@echo "Singularity image: $(IMAGE)_$(TAG).sif"

clean:
	docker rmi $(IMAGE):$(TAG) 2>/dev/null || true
