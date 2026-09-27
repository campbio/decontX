# Standard targets expected by r-bioc-dev-standards, plus decontX's
# coverage/build/clean. Run `make help` to list them. Recipe lines must
# start with a tab. All Suggests must be installed before running checks.

.PHONY: help docs test test-one check check-full bioccheck lint coverage coverage-report build clean

help:  ## List targets
	@grep -E '^[a-zA-Z_-]+:.*## ' $(MAKEFILE_LIST) | \
	  awk 'BEGIN {FS = ":.*## "}; {printf "  %-12s %s\n", $$1, $$2}'

docs:  ## Regenerate man/*.Rd and NAMESPACE from roxygen comments
	Rscript -e 'devtools::document()'

test:  ## Run the full test suite
	Rscript -e 'devtools::test(stop_on_failure = TRUE)'

test-one:  ## Run matching test files: make test-one FILTER=<pattern>
	@test -n "$(FILTER)" || { echo "Usage: make test-one FILTER=<pattern>"; exit 1; }
	Rscript -e 'devtools::test(filter = "$(FILTER)", stop_on_failure = TRUE)'

check:  ## Quick R CMD check: skips vignettes and the PDF manual
	Rscript -e 'rcmdcheck::rcmdcheck(args = c("--no-manual", "--ignore-vignettes"), build_args = "--no-build-vignettes", error_on = "warning", check_dir = tempdir())'

check-full:  ## Full check: rebuilds vignettes, runs \donttest examples; needs all Suggests
	Rscript -e 'devtools::check(document = FALSE, vignettes = TRUE, run_dont_test = TRUE, force_suggests = TRUE, error_on = "warning", check_dir = tempdir())'

bioccheck:  ## BiocCheckGitClone on the repo, then BiocCheck on a built tarball
	Rscript -e 'BiocCheck::BiocCheckGitClone(".", "quit-with-status" = TRUE)'
	@tmp="$$(mktemp -d)"; \
	  cd "$$tmp" && R CMD build "$(CURDIR)" && \
	  Rscript -e 'tb <- list.files(pattern = "[.]tar[.]gz$$")[1]; message(sprintf("Tarball %s: %.2f MB (limit 10 MB)", tb, file.size(tb) / 1e6)); BiocCheck::BiocCheck(tb, "quit-with-status" = TRUE)'; \
	  status=$$?; echo "BiocCheck output kept in $$tmp"; exit $$status

lint:  ## Run lintr on the package (reports only; changes nothing)
	Rscript -e 'print(lintr::lint_package())'

coverage:  ## Print test coverage, overall and per file
	Rscript -e 'print(covr::package_coverage())'

coverage-report:  ## HTML coverage report (opens a browser)
	Rscript -e 'covr::report(covr::package_coverage())'

build:  ## Build the source tarball in the repo root
	R CMD build .

clean:  ## Remove build artifacts (compiled objects, tarballs)
	rm -f src/*.o src/*.so src/*.dll decontX_*.tar.gz
