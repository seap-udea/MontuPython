.PHONY: help show status clean cleanall cleancrap cleanout cleandist cleanenv cleangit cleandocs \
	install install-dev build env addall commit pull push release notebooks readme import \
	docs docs-install docs-prepare docs-build docs-clean \
	test-install test test-docstrings test-notebooks test-structure test-desktop test-organic \
	desktop-show desktop-install-build desktop-build desktop-clean \
	desktop-package desktop-release desktop-ci subsets

##################################################################
# VARIABLES
##################################################################
SHELL := /bin/bash
BRANCH := $(shell git rev-parse --abbrev-ref HEAD)
VERSION := $(shell tail -n 1 .versions)
COMMIT_MSG ?= [MAN] Maintenance
RELMODE ?= release
PYTHON ?= python3
PIP ?= pip3
PACKNAME := montu
SPHINXOPTS ?=
DESKTOP_VERSION := $(shell grep '^version' montu_gui/version.py | head -1 | sed "s/.*['\"]\\(.*\\)['\"].*/\\1/")

help:
	@echo "MontuPython Development Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  show        - Show current version and branch"
	@echo "  status      - Show git status"
	@echo "  install     - Install the package"
	@echo "  install-dev - Install editable package + dev, test, and desktop build deps"
	@echo "  build       - Build distribution packages"
	@echo "  env         - Create local dev environment (.venv)"
	@echo "  clean       - Remove cache files"
	@echo "  cleanall    - Deep clean (build + caches)"
	@echo "  push        - Commit and push current branch"
	@echo "  release     - Release a new version (make release RELMODE=release VERSION=x.y.z)"
	@echo "  release-pipeline - Full release workflow (./release-pipeline.sh --version x.y.z [--tag a.b.c])"
	@echo "  notebooks     - Execute README + examples notebooks (embed plots)"
	@echo "  readme      - Convert README.ipynb to README.md"
	@echo "  docs         - Prepare examples and build Sphinx HTML (prepare + build)"
	@echo "  docs-install - Install documentation dependencies (sphinx, myst-parser, nbsphinx, etc.)"
	@echo "  docs-prepare - Copy example notebooks into docs/examples/ and update examples.rst"
	@echo "  docs-build   - Prepare examples and build Sphinx HTML documentation"
	@echo "  docs-clean   - Remove the docs/_build directory"
	@echo "  test-install - Install test dependencies"
	@echo "  test         - Run the full pytest suite"
	@echo "  test-docstrings - Run tests derived from docstring examples"
	@echo "  test-notebooks  - Run tests derived from example notebooks"
	@echo "  test-structure  - Validate example notebook structure"
	@echo "  test-desktop    - Run MontuPython Desktop module and example tests"
	@echo "  test-organic    - Run organic ephemeris/stellar regression snapshots"
	@echo "  organic-snapshots - Regenerate organic CSV regression references"
	@echo "  pin-astronomy   - Pin ephem/pymeeus/pyplanets and refresh organic snapshots (VERSION=x.y.z)"
	@echo ""
	@echo "MontuPython Desktop ($(DESKTOP_VERSION)):"
	@echo "  desktop-show          - Show desktop app version"
	@echo "  desktop-install-build - Install PyInstaller build deps in .desktop-build"
	@echo "  desktop-build         - Build .app (macOS) or onedir bundle"
	@echo "  desktop-clean         - Remove desktop build artifacts"
	@echo "  desktop-package       - Build and create zip/dmg (macOS) or zip (Windows)"
	@echo "  desktop-release       - Bump version + build (make desktop-release VERSION=0.1.2)"
	@echo "  desktop-ci            - Trigger GitHub CI build (make desktop-ci TAG=desktop-v0.1.2)"

show:
	@echo "Version: $(VERSION)"
	@echo "Branch: $(BRANCH)"

status:
	@git status

##################################################################
# BASIC RULES
##################################################################
clean: cleancrap

cleanall: _dev_cleanall cleancrap cleanout cleandist cleanenv cleangit cleandocs

#=========================
# Clean
#=========================
cleancrap:
	@echo "Cleaning crap..."
	@-find . -name "*~" -delete
	@-find . -name "#*#" -delete
	@-find . -name "#*" -delete
	@-find . -name ".#*" -delete
	@-find . -name ".#*#" -delete
	@-find . -name ".DS_Store" -delete
	@-find . -name "Icon*" -delete
	@-find . -name "montu_dem*" -type d | xargs rm -fr
	@-find . -name "*.egg-info*" -type d | xargs rm -fr

cleanout:
	@echo "Cleaning all compiled objects..."
	@-find . -name "*.o" -not -path "./.venv/*" -not -path "./.desktop-build/*" -delete
	@-find . -name "*.opp" -not -path "./.venv/*" -not -path "./.desktop-build/*" -delete
	@-find . -name "*.gcno" -not -path "./.venv/*" -not -path "./.desktop-build/*" -delete
	@-find . -name "*.gcda" -not -path "./.venv/*" -not -path "./.desktop-build/*" -delete
	@-find . -name "*.gcov" -not -path "./.venv/*" -not -path "./.desktop-build/*" -delete
	@-find . -name "*.info" -not -path "./.venv/*" -not -path "./.desktop-build/*" -delete
	@-find . -name "*.out" -not -path "./.venv/*" -not -path "./.desktop-build/*" -delete
	@-find . -name "*.tout" -not -path "./.venv/*" -not -path "./.desktop-build/*" -delete
	@-find . -name "*.so" -not -path "./.venv/*" -not -path "./.desktop-build/*" -delete
	@-find . -name ".ipynb_checkpoints" -not -path "./.venv/*" -not -path "./.desktop-build/*" -type d | xargs rm -fr
	@-find . -name "__pycache__" -not -path "./.venv/*" -not -path "./.desktop-build/*" -type d | xargs rm -fr

cleandist:
	@-rm -rf dist/
	@-rm -rf build/
	@-rm -rf "dist/MontuPython Desktop.app" dist/MontuPython-Desktop dist/desktop

cleanenv:
	@echo "Cleaning local environments..."
	@-rm -rf .venv
	@-rm -rf .desktop-build

cleandocs:
	@echo "Cleaning documentation build..."
	@-rm -rf docs/_build

##################################################################
# PACKAGE RULES
##################################################################
install:
	$(PYTHON) -m pip install .

install-dev: env
	@echo "Installing development requirements in .venv..."
	@if [ -f requirements.txt ]; then . .venv/bin/activate && pip install -r requirements.txt; fi
	@if [ -f requirements-test.txt ]; then . .venv/bin/activate && pip install -r requirements-test.txt; fi
	@if [ -f requirements-desktop-build.txt ]; then . .venv/bin/activate && pip install -r requirements-desktop-build.txt; fi

env:
	@echo "Creating local development environment..."
	@test -d .venv || $(PYTHON) -m venv .venv
	@echo "Installing dependencies from setup.py..."
	@. .venv/bin/activate && pip install --upgrade pip
	@. .venv/bin/activate && pip install -e .
	@echo "______________________________________________________________________"
	@echo "Environment setup complete."
	@echo "To activate the environment, run:"
	@echo "source .venv/bin/activate"

build: clean
	$(PYTHON) -m build

##################################################################
# GIT
##################################################################
cleangit:
	@echo "Purging old and heavy files from local git history..."
	@git reflog expire --expire=now --all
	@git gc --prune=now --aggressive

addall: cleanall
	@echo "Adding..."
	@-git add -A .

commit:
	@echo "Committing..."
	@git commit -am "$(COMMIT_MSG)"
	@-git push origin $(BRANCH)

pull:
	@echo "Pulling new files..."
	@-git pull origin $(BRANCH)

push:
	@echo "Committing tracked changes (if any)..."
	@if ! git diff --quiet || ! git diff --cached --quiet || [ -n "$$(git status --porcelain)" ]; then \
		git add . && \
		files="$$(git diff --cached --name-only | paste -sd', ' - || true)" && \
		msg="$(COMMIT_MSG)" && \
		if [ "$(origin COMMIT_MSG)" != "command line" ] && [ "$(origin COMMIT_MSG)" != "environment" ]; then \
			if [ -n "$$files" ]; then msg="$$msg [$$files]"; fi; \
		fi && \
		git commit -m "$$msg"; \
	else \
		echo "Working tree is clean (tracked files); nothing to commit."; \
	fi
	@echo "Pushing current branch..."
	@git push -u origin HEAD

##################################################################
# RELEASE
##################################################################
# Example: make release RELMODE=release VERSION=0.9.11
release:
	@echo "Releasing a new version..."
	@bash bin/release.sh $(RELMODE) $(VERSION)

notebooks:
	@echo "Executing README and example notebooks..."
	@$(PYTHON) bin/execute_notebooks.py

readme:
	$(PYTHON) -m nbconvert README.ipynb --to markdown

import:
	@$(PYTHON) -c "from montu import *;print(version)"

##################################################################
# DOCUMENTATION
##################################################################
docs: docs-prepare docs-build

docs-prepare:
	@echo "Preparing example notebooks..."
	@$(PYTHON) bin/examples_doc.py
	@echo "Injecting Plotly HTML for Sphinx compatibility..."
	@$(PYTHON) bin/fix_notebook_plotly.py docs/examples/*.ipynb

docs-install:
	@echo "Installing documentation dependencies..."
	@$(PYTHON) -m pip install -r docs/requirements.txt

docs-build: docs-prepare
	@echo "Building Sphinx HTML documentation..."
	@cd docs && $(PYTHON) -m sphinx -M html . _build $(SPHINXOPTS)
	@echo "Documentation available at docs/_build/html/index.html"

docs-clean:
	@echo "Cleaning documentation build..."
	@-rm -rf docs/_build

##################################################################
# TESTS
##################################################################
test-install:
	@echo "Installing test dependencies..."
	@$(PYTHON) -m pip install -r requirements-test.txt

test:
	@echo "Running full test suite..."
	@$(PYTHON) -m pytest

test-docstrings:
	@echo "Running tests derived from docstring examples..."
	@$(PYTHON) -m pytest -m docstrings

test-notebooks:
	@echo "Running tests derived from example notebooks..."
	@$(PYTHON) -m pytest -m notebooks

test-structure:
	@echo "Validating example notebook structure..."
	@$(PYTHON) -m pytest -m structure

test-desktop:
	@echo "Running MontuPython Desktop tests..."
	@$(PYTHON) -m pytest montu_gui/tests -m desktop

test-organic:
	@echo "Running organic snapshot regression tests..."
	@$(PYTHON) -m pytest montu/tests/test_planetary_ephemeris_organic.py montu/tests/test_stellar_positions_organic.py

organic-snapshots:
	@echo "Regenerating organic ephemeris and stellar position snapshots..."
	@$(if $(wildcard .venv/bin/python),.venv/bin/python,$(PYTHON)) bin/generate_organic_stellar_positions.py
	@$(if $(wildcard .venv/bin/python),.venv/bin/python,$(PYTHON)) bin/generate_organic_planetary_ephemeris.py

pin-astronomy:
	@if [ -z "$(VERSION)" ]; then \
		echo "Usage: make pin-astronomy VERSION=x.y.z"; \
		exit 1; \
	fi
	@echo "Pinning astronomy stack for release $(VERSION)..."
	@$(if $(wildcard .venv/bin/python),.venv/bin/python,$(PYTHON)) bin/update_release_astronomy_stack.py --version "$(VERSION)" --python $(if $(wildcard .venv/bin/python),.venv/bin/python,$(PYTHON))

##################################################################
# MONTUPYTHON DESKTOP
##################################################################
desktop-show:
	@echo "Desktop version: $(DESKTOP_VERSION)"

desktop-install-build:
	@test -d .desktop-build || $(PYTHON) -m venv .desktop-build
	@. .desktop-build/bin/activate && pip install -q --upgrade pip
	@. .desktop-build/bin/activate && pip install -q -e .
	@. .desktop-build/bin/activate && pip install -q -r requirements.txt
	@. .desktop-build/bin/activate && pip install -q -r requirements-desktop-build.txt
	@echo "Desktop build environment ready (.desktop-build)."

desktop-build:
	@bash bin/build-desktop.sh --no-package

desktop-clean:
	@echo "Cleaning desktop build artifacts..."
	@-rm -rf dist/desktop dist/MontuPython-Desktop "dist/MontuPython Desktop.app" build/montu-desktop

desktop-package:
	@bash bin/build-desktop.sh

desktop-release:
	@if [ -z "$(VERSION)" ]; then \
		echo "Usage: make desktop-release VERSION=x.y.z"; \
		exit 1; \
	fi
	@bash bin/desktop-release.sh "$(VERSION)"

# Push a desktop tag to trigger GitHub Actions (Mac + Windows builds).
# Example: make desktop-ci TAG=desktop-v0.1.2
desktop-ci:
	@if [ -z "$(TAG)" ]; then \
		echo "Usage: make desktop-ci TAG=desktop-v0.1.2"; \
		exit 1; \
	fi
	@git push origin "$(TAG)"

subsets:
	@echo "Generating stellar subsets..."
	$(PYTHON) bin/cat_subsets.py
# --- dev/cleanall (auto) ---
include .dev_common.mk
