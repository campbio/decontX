# Makefile for an r-bioc-dev-standards package. Run `make help` to list
# targets. Recipe lines must start with a tab.
#
# The standard targets (test, check, bioccheck, ...) live in the shared
# standards.mk in r-bioc-dev-standards, not here. It's downloaded to a cache
# on first use and refreshed at every Claude session start.

# Package settings. Uncomment to change a default (see standards.mk).
# FORCE_SUGGESTS = FALSE

R_BIOC_STANDARDS_REF ?= v1
R_BIOC_STANDARDS_BASE ?= https://raw.githubusercontent.com/campbio/r-bioc-dev-standards/$(R_BIOC_STANDARDS_REF)
STANDARDS_MK := $(or $(XDG_CACHE_HOME),$(HOME)/.cache)/r-bioc-dev-standards/$(R_BIOC_STANDARDS_REF)/standards.mk

include $(STANDARDS_MK)

$(STANDARDS_MK):
	@mkdir -p "$(@D)"
	curl -fsSL --max-time 30 "$(R_BIOC_STANDARDS_BASE)/shared/standards.mk" -o "$@.tmp"
	@mv -f "$@.tmp" "$@"

# Extra targets for this package go below, each listed in AGENTS.md. To
# replace a standard target, define it here; make warns that it overrides
# the shared recipe.

.PHONY: coverage-report build clean

coverage-report:  ## HTML coverage report (opens a browser)
	Rscript -e 'covr::report(covr::package_coverage())'

build:  ## Build the source tarball in the repo root
	R CMD build .

clean:  ## Remove build artifacts (compiled objects, tarballs)
	rm -f src/*.o src/*.so src/*.dll decontX_*.tar.gz
