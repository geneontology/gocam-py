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

# note: "help" MUST be the first target in the file,
# when the user types "make" they should get help info
.PHONY: help
help:
	@echo ""
	@echo "make install -- install dependencies"
	@echo "make audit -- scan dependencies for malware and known vulnerabilities"
	@echo "make gen-project -- generate LinkML project artifacts"
	@echo "make gendoc -- generate schema documentation"
	@echo "make gen-examples -- refresh committed model examples"
	@echo "make gen-test-inputs -- refresh committed test inputs"
	@echo "make test -- run all tests"
	@echo "make test-python -- run Python tests"
	@echo "make lint -- lint the LinkML schema"
	@echo "make lint-python -- check Python formatting and linting"
	@echo "make lint-fix-python -- fix Python formatting and linting"
	@echo "make type-check -- type-check Python code"
	@echo "make site -- generate the local documentation site"
	@echo "make serve -- serve documentation locally"
	@echo "make testdoc -- generate and serve documentation locally"
	@echo "make deploy -- deploy the documentation site"
	@echo "make clean -- remove generated artifacts"
	@echo "make help -- show this help"
	@echo ""

# install any dependencies required for building
.PHONY: install
install:
	UV_MALWARE_CHECK=1 uv sync --frozen --all-extras

.PHONY: audit
audit: install
	uv audit --frozen

.PHONY: all
all: site

.PHONY: site
site: gen-project gendoc $(PYDANTIC)

.PHONY: deploy
deploy: all mkd-gh-deploy

.PHONY: gen-examples
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

.PHONY: gen-project
gen-project:
	$(RUN) gen-project ${CONFIG_YAML} --exclude excel --exclude graphql -d $(DEST) $(SOURCE_SCHEMA_PATH)

.PHONY: test
test: test-schema test-python test-examples

.PHONY: test-schema
test-schema:
	$(RUN) gen-project ${CONFIG_YAML} --exclude excel --exclude graphql -d tmp $(SOURCE_SCHEMA_PATH)

.PHONY: test-python
test-python:
	$(RUN) pytest tests

.PHONY: lint
lint:
	$(RUN) linkml-lint $(SOURCE_SCHEMA_PATH)

.PHONY: lint-python
lint-python:
	$(RUN) ruff format --check
	$(RUN) ruff check

.PHONY: lint-fix-python
lint-fix-python:
	$(RUN) ruff check --fix
	$(RUN) ruff format

.PHONY: type-check
type-check:
	$(RUN) ty check

.PHONY: test-examples
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
.PHONY: serve
serve: mkd-serve

# Python datamodel
$(PYMODEL):
	mkdir -p $@

$(PYDANTIC): $(SOURCE_SCHEMA_PATH)
	$(RUN) linkml generate pydantic $< > $@.tmp && mv $@.tmp $@

$(DOCDIR):
	mkdir -p $@

.PHONY: gendoc
gendoc: $(DOCDIR)
	cp -rf src/docs/* $(DOCDIR) ; \
	$(RUN) gen-doc -d $(DOCDIR) $(SOURCE_SCHEMA_PATH)

.PHONY: testdoc
testdoc: gendoc serve

MKDOCS = $(RUN) mkdocs
mkd-%:
	$(MKDOCS) $*

.PHONY: clean
clean:
	rm -rf $(DEST)
	rm -rf tmp
	rm -fr docs/*
	rm -fr $(PYDANTIC)
