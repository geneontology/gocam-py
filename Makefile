MAKEFLAGS += --warn-undefined-variables
SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DEFAULT_GOAL := help
.DELETE_ON_ERROR:
.SUFFIXES:
.SECONDARY:

RUN := uv run
SOURCE_SCHEMA_PATH := src/gocam/schema/gocam.yaml
DEST := project
PYMODEL := src/gocam/datamodel
PYDANTIC := $(PYMODEL)/gocam.py
DOCDIR := docs
CONFIG_YAML := --config-file config.yaml

.PHONY: all clean gen-project gen-examples gendoc lint lint-python lint-fix-python audit install

# note: "help" MUST be the first target in the file,
# when the user types "make" they should get help info
help:
	@echo ""
	@echo "make site -- makes site locally"
	@echo "make install -- install dependencies"
	@echo "make audit -- scan dependencies for malware and known vulnerabilities"
	@echo "make test -- runs tests"
	@echo "make lint -- perform linting"
	@echo "make testdoc -- builds docs and runs local test server"
	@echo "make deploy -- deploys site"
	@echo "make translate-collection -- translate GO-CAM collection to networkx and cx2"
	@echo "make help -- show this help"
	@echo ""

# install any dependencies required for building
install:
	UV_MALWARE_CHECK=1 uv sync --frozen --all-extras
.PHONY: install

audit: install
	uv audit --frozen

# EXPERIMENTAL
create-data-harmonizer:
	npm init data-harmonizer $(SOURCE_SCHEMA_PATH)

all: site
site: gen-project gendoc $(PYDANTIC)
%.yaml: gen-project
deploy: all mkd-gh-deploy

compile-sheets:
	$(RUN) sheets2linkml --gsheet-id $(SHEET_ID) $(SHEET_TABS) > $(SHEET_MODULE_PATH).tmp && mv $(SHEET_MODULE_PATH).tmp $(SHEET_MODULE_PATH)

gen-examples:
	$(RUN) gocam fetch --format yaml 663d668500002178 > src/data/examples/Model-663d668500002178.yaml
	$(RUN) gocam fetch --format json 663d668500002178 > src/data/examples/Model-663d668500002178.json

.PHONY: gen-test-inputs
gen-test-inputs:
	$(RUN) gocam fetch --format yaml 63f809ec00000701 > tests/input/Model-63f809ec00000701.yaml
	$(RUN) gocam fetch --format yaml 568b0f9600000284 > tests/input/Model-568b0f9600000284.yaml
	$(RUN) gocam fetch --format yaml 663d668500002178 > tests/input/Model-663d668500002178.yaml
	$(RUN) gocam fetch --format yaml 6606056e00002011 > tests/input/Model-6606056e00002011.yaml

# generates all project files

gen-project:
	$(RUN) gen-project ${CONFIG_YAML} --exclude excel --exclude graphql -d $(DEST) $(SOURCE_SCHEMA_PATH)

test: test-schema test-python test-examples

test-schema:
	$(RUN) gen-project ${CONFIG_YAML} --exclude excel --exclude graphql -d tmp $(SOURCE_SCHEMA_PATH)

test-python:
	$(RUN) pytest tests

lint:
	$(RUN) linkml-lint $(SOURCE_SCHEMA_PATH)

lint-python:
	$(RUN) ruff format --check
	$(RUN) ruff check

lint-fix-python:
	$(RUN) ruff check --fix
	$(RUN) ruff format

type-check:
	$(RUN) ty check

convert-examples-to-%:
	$(patsubst %, $(RUN) linkml-convert  % -s $(SOURCE_SCHEMA_PATH) -C Person, $(shell ${SHELL} find src/data/examples -name "*.yaml"))

examples/%.yaml: src/data/examples/%.yaml
	$(RUN) linkml-convert -s $(SOURCE_SCHEMA_PATH) -C Person $< -o $@
examples/%.json: src/data/examples/%.yaml
	$(RUN) linkml-convert -s $(SOURCE_SCHEMA_PATH) -C Person $< -o $@
examples/%.ttl: src/data/examples/%.yaml
	$(RUN) linkml-convert -P EXAMPLE=http://example.org/ -s $(SOURCE_SCHEMA_PATH) -C Person $< -o $@

test-examples: examples/output

examples/output: src/gocam/schema/gocam.yaml
	mkdir -p $@
	$(RUN) linkml-run-examples \
		--output-formats json \
		--output-formats yaml \
		--counter-example-input-directory src/data/examples/invalid \
		--input-directory src/data/examples/valid \
		--output-directory $@ \
		--schema $< > $@/README.md

# Test documentation locally
serve: mkd-serve

# Python datamodel
$(PYMODEL):
	mkdir -p $@

$(PYDANTIC): $(SOURCE_SCHEMA_PATH)
	$(RUN) linkml generate pydantic $< > $@.tmp && mv $@.tmp $@

$(DOCDIR):
	mkdir -p $@

gendoc: $(DOCDIR)
	cp -rf src/docs/* $(DOCDIR) ; \
	$(RUN) gen-doc -d $(DOCDIR) $(SOURCE_SCHEMA_PATH)

testdoc: gendoc serve

MKDOCS = $(RUN) mkdocs
mkd-%:
	$(MKDOCS) $*

# Translate GO-CAM collection to networkx and cx2 formats
translate-collection:
	$(RUN) gocam translate-collection

# Translate GO-CAM collection with custom parameters (example)
translate-collection-test:
	$(RUN) gocam -v translate-collection --limit 5

clean:
	rm -rf $(DEST)
	rm -rf tmp
	rm -fr docs/*
	rm -fr $(PYDANTIC)
