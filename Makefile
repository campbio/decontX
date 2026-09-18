# Canonical developer commands for decontX — one answer to "how do I test this"
# for humans, agents, and CI. Non-obvious flags live here, not in prose.
#
# decontX has no Shiny app and no pkgdown site, so there are no app/site
# targets (see AGENTS.md).

.PHONY: test check bioccheck docs lint build clean

test:        ## fast loop — run after every change
	Rscript -e 'devtools::test()'

check:       ## full check — run before opening a PR
	_R_CHECK_FORCE_SUGGESTS_=false Rscript -e 'rcmdcheck::rcmdcheck(args = c("--no-manual"), error_on = "warning")'

bioccheck:   ## build tarball first — matches how the Bioc build system runs it
	R CMD build . && Rscript -e 'BiocCheck::BiocCheck(Sys.glob("decontX_*.tar.gz")); BiocCheck::BiocCheckGitClone(".")'

docs:        ## regenerate man/ and NAMESPACE from roxygen comments
	Rscript -e 'devtools::document()'

lint:
	Rscript -e 'lintr::lint_package()'

build:
	R CMD build .

clean:       ## remove build artifacts (compiled objects, tarballs)
	rm -f src/*.o src/*.so src/*.dll decontX_*.tar.gz
