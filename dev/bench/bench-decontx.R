# Benchmark speed and memory of decontX() and decontPro().
#
# Purpose
#   Times decontX() (single batch with user or estimated clusters, four
#   batches, and with a background matrix) and decontPro() on simulated
#   data, profiles one mid-size run of each, and benchmarks real PBMC
#   data. Used for the speed and memory work (phase I-e of the deep
#   review plan): each optimization PR reports before/after numbers from
#   this script.
#
# How to run (from the repository root)
#   Rscript dev/bench/bench-decontx.R [mode] [outdir]
#
#   mode    "sim" (default): simulated data at several sizes.
#           "real": 10x PBMC data from TENxPBMCData (downloads about
#           30 MB into the ExperimentHub cache on first use).
#           "profile": Rprof (and profvis, if installed) on one mid-size
#           decontX() run with estimated clusters, with user clusters,
#           and with four batches, and one decontPro() run.
#   outdir  Where results go. Default: a new folder in tempdir(), which
#           is deleted when the R session ends, so pass a folder to keep
#           the results.
#
#   Keep the machine awake while it runs (on macOS:
#   caffeinate -dims Rscript dev/bench/bench-decontx.R ...); otherwise
#   wall times include time asleep. A median_s much larger than cpu_s in
#   the results means the run was suspended or starved: rerun it. Run
#   nothing else heavy at the same time, and do not edit this script
#   while it runs (R reads it lazily).
#   With the defaults, sim mode takes about an hour; decontPro alone is
#   several minutes per size.
#
# Installed build only
#   devtools::load_all() compiles the C++ code without optimization and
#   runs the EM step about 9 times slower, so this script never uses it.
#   It copies the working tree (tracked files plus untracked files that
#   are not ignored) to a temporary folder, installs it into
#   <outdir>/lib with R CMD INSTALL, and loads decontX from there in a
#   fresh R subprocess for every case. The first install compiles the
#   Stan model and takes several minutes.
#
# Environment variables
#   DECONTX_BENCH_LIB        Library that already holds an installed
#                            decontX. Skips the install; also use it to
#                            point at another commit's install for
#                            before/after comparisons.
#   DECONTX_BENCH_LABEL      Label stored in the results (default: the
#                            short git commit of the working tree).
#   DECONTX_BENCH_SIZES      Cell counts for sim mode
#                            (default "10000,50000,100000").
#   DECONTX_BENCH_RUNS       decontX runs for sim mode (default
#                            "z_supplied,z_estimated,batched4,background").
#   DECONTX_BENCH_PRO_CELLS  Cell counts for decontPro in sim mode
#                            (default "2000,5000"; "" skips decontPro).
#                            Profile mode uses the largest.
#   DECONTX_BENCH_ITER       Iterations per case (default: 3 up to 10k
#                            cells, 2 up to 50k, then 1).
#   DECONTX_BENCH_PROFILE_CELLS  Cells for profile mode (default 50000).
#   DECONTX_BENCH_PBMC       TENxPBMCData data set for real mode
#                            (default "pbmc4k").
#
# Outputs (in outdir)
#   bench-<mode>-<time>.csv   One row per case: wall time (median, min,
#                             max), CPU time, R-allocated memory, peak
#                             RSS of the case's process, EM iterations.
#   sessionInfo-<mode>.txt    sessionInfo() of a worker with decontX
#                             loaded from the benchmark library.
#   profile/                  Rprof output, by.self/by.total tables, a
#                             per-phase table, and profvis HTML.
#   data/                     Cached simulated data (dgCMatrix, RDS).
#
# Memory
#   bench::bench_memory() counts memory allocated by R only. Memory
#   allocated by C++ (Eigen, Stan, scrapper) is not counted, so each case
#   also runs in its own process under /usr/bin/time, which reports the
#   peak resident set size (RSS) of that process. Peak RSS includes R
#   itself and the loaded input data; rss_before_mb is the RSS just
#   before the timed call.
#
#   simulateContamination() builds dense integer genes x cells matrices
#   (about 0.8 GB each at 2,000 genes x 100,000 cells) and peaks at a
#   few GB, so data are simulated once per size in a separate process
#   and cached as a dgCMatrix.
#
# Dev-only tools
#   bench (required) and profvis (optional, profile mode) are
#   development tools, not DESCRIPTION dependencies. Install them in
#   your own library. TENxPBMCData is needed for real mode only.


## ---------------------------------------------------------------------
## Settings
## ---------------------------------------------------------------------

## Simulation settings shared by all sizes
simGenes <- 2000
simClusters <- 10
simCountRange <- c(1000, 4000)
simSeed <- 12345
simBatches <- 4

## decontPro panel: an ADT panel of this many antibodies
proAdts <- 40
proCellTypes <- 5

## Read an environment variable, with a default when it is unset. A
## variable set to "" stays "" (so DECONTX_BENCH_PRO_CELLS="" skips
## decontPro).
envOr <- function(name, default) {
  value <- Sys.getenv(name, unset = NA)
  if (is.na(value)) default else value
}

## Split "a,b,c" into a character vector, dropping empty entries
splitList <- function(x) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  parts[nzchar(parts)]
}

## Path of this script, so the parent can start workers with it
scriptPath <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  fileArg <- grep("^--file=", args, value = TRUE)
  if (length(fileArg) == 0) {
    stop("Run this script with Rscript.")
  }
  normalizePath(sub("^--file=", "", fileArg[1]))
}

## Iterations for a case with `cells` cells
iterationsFor <- function(cells) {
  fixed <- Sys.getenv("DECONTX_BENCH_ITER", unset = "")
  if (nzchar(fixed)) {
    return(as.integer(fixed))
  }
  if (cells <= 10000) {
    3L
  } else if (cells <= 50000) {
    2L
  } else {
    1L
  }
}


## ---------------------------------------------------------------------
## Installing the working tree into a temporary library
## ---------------------------------------------------------------------

## Copy the files git knows about (tracked, plus untracked files that are
## not ignored) so the install sees uncommitted edits but no build
## leftovers, then install the copy. Installing the repository folder
## directly would leave compiled objects in src/.
installWorkingTree <- function(repoDir, lib) {
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  files <- system2("git",
    c("-C", shQuote(repoDir), "ls-files", "--cached", "--others",
      "--exclude-standard"),
    stdout = TRUE
  )
  files <- files[file.exists(file.path(repoDir, files))]
  ## Skip dev/ and other build-ignored folders the install does not need
  files <- files[!grepl("^(dev|\\.github|\\.claude)/", files)]

  copyDir <- file.path(tempfile("decontx-src-"), "decontX")
  for (subDir in unique(dirname(files))) {
    dir.create(file.path(copyDir, subDir),
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  copied <- file.copy(file.path(repoDir, files), file.path(copyDir, files))
  if (!all(copied)) {
    stop("Could not copy the working tree to ", copyDir)
  }

  ## Parallel make speeds up compiling the Stan model
  if (!nzchar(Sys.getenv("MAKEFLAGS"))) {
    Sys.setenv(MAKEFLAGS = "-j4")
  }
  message("Installing decontX from the working tree into ", lib)
  status <- system2(file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-docs", "--no-multiarch", "--no-test-load",
      paste0("--library=", shQuote(lib)), shQuote(copyDir)),
    stdout = file.path(lib, "install.log"),
    stderr = file.path(lib, "install.log")
  )
  if (status != 0) {
    stop("R CMD INSTALL failed; see ", file.path(lib, "install.log"))
  }
  unlink(dirname(copyDir), recursive = TRUE)
  invisible(lib)
}

## Return the benchmark library, installing decontX into it if needed
benchLibrary <- function(repoDir, outdir) {
  lib <- Sys.getenv("DECONTX_BENCH_LIB", unset = "")
  if (nzchar(lib)) {
    if (!file.exists(file.path(lib, "decontX", "DESCRIPTION"))) {
      stop("DECONTX_BENCH_LIB has no installed decontX: ", lib)
    }
    message("Using the decontX already installed in ", lib)
    return(normalizePath(lib))
  }
  lib <- file.path(outdir, "lib")
  installWorkingTree(repoDir, lib)
  normalizePath(lib)
}


## ---------------------------------------------------------------------
## Running a case in its own process
## ---------------------------------------------------------------------

## Peak RSS in MB from the output of /usr/bin/time (macOS -l prints
## bytes, GNU time -v prints kilobytes)
parsePeakRss <- function(lines) {
  macLine <- grep("maximum resident set size", lines, value = TRUE)
  if (length(macLine) > 0) {
    return(as.numeric(sub("^\\s*(\\d+).*$", "\\1", macLine[1])) / 2^20)
  }
  gnuLine <- grep("Maximum resident set size", lines, value = TRUE)
  if (length(gnuLine) > 0) {
    return(as.numeric(sub("^.*:\\s*(\\d+).*$", "\\1", gnuLine[1])) / 2^10)
  }
  NA_real_
}

## Start `Rscript <this script> worker <task> <outdir> <lib>` and return
## the worker's result (a one-row data frame) with its peak RSS added
runWorker <- function(task, outdir, lib) {
  resultFile <- file.path(outdir, "work", paste0(task, ".rds"))
  logFile <- file.path(outdir, "work", paste0(task, ".log"))
  dir.create(dirname(resultFile), showWarnings = FALSE, recursive = TRUE)
  unlink(resultFile)

  rscript <- file.path(R.home("bin"), "Rscript")
  args <- c(
    shQuote(scriptPath()), "worker", task, shQuote(outdir), shQuote(lib)
  )
  timeFlag <- if (Sys.info()[["sysname"]] == "Darwin") "-l" else "-v"
  useTime <- file.exists("/usr/bin/time")
  ## Put the benchmark library first so decontX loads from it
  libPaths <- paste(c(lib, .libPaths()), collapse = .Platform$path.sep)

  message(format(Sys.time(), "%H:%M:%S"), "  running ", task)
  if (useTime) {
    status <- system2("/usr/bin/time", c(timeFlag, rscript, args),
      stdout = logFile, stderr = logFile,
      env = paste0("R_LIBS=", shQuote(libPaths))
    )
  } else {
    status <- system2(rscript, args,
      stdout = logFile, stderr = logFile,
      env = paste0("R_LIBS=", shQuote(libPaths))
    )
  }
  if (status != 0 || !file.exists(resultFile)) {
    warning("Task ", task, " failed; see ", logFile, call. = FALSE)
    return(NULL)
  }
  result <- readRDS(resultFile)
  if (is.data.frame(result)) {
    result$peak_rss_mb <- parsePeakRss(readLines(logFile))
  }
  result
}


## ---------------------------------------------------------------------
## Data
## ---------------------------------------------------------------------

## Simulate `cells` cells with simulateContamination() and save counts
## (dgCMatrix), true clusters, batch labels and a background matrix
simulateData <- function(cells, file) {
  sim <- decontX::simulateContamination(
    C = cells,
    G = simGenes,
    K = simClusters,
    NRange = simCountRange,
    seed = simSeed
  )
  counts <- methods::as(sim$observedCounts, "CsparseMatrix")
  counts <- methods::as(counts, "dMatrix")
  z <- sim$z
  rm(sim)
  invisible(gc())

  ## Contiguous batches of equal size, as from separate samples
  batchSize <- ceiling(cells / simBatches)
  batchIndex <- rep(seq_len(simBatches), each = batchSize)[seq_len(cells)]
  batch <- paste0("batch", batchIndex)

  ## Background (empty) droplets: as many as cells, 10-200 counts each,
  ## drawn from the overall gene profile of the cells
  background <- withr::with_seed(simSeed, {
    profile <- Matrix::rowSums(counts) / sum(counts)
    dropletCounts <- sample(10:200, cells, replace = TRUE)
    genes <- sample.int(nrow(counts), sum(dropletCounts),
      replace = TRUE, prob = profile
    )
    Matrix::sparseMatrix(
      i = genes,
      j = rep(seq_len(cells), dropletCounts),
      x = 1,
      dims = c(nrow(counts), cells),
      dimnames = list(rownames(counts), paste0("Empty_", seq_len(cells)))
    )
  })

  saveRDS(
    list(counts = counts, z = z, batch = batch, background = background),
    file
  )
}

## Simulate an ADT panel for decontPro: few features, dense counts
simulateAdt <- function(cells) {
  sim <- decontX::simulateContamination(
    C = cells,
    G = proAdts,
    K = proCellTypes,
    NRange = c(1000, 5000),
    seed = simSeed
  )
  list(counts = sim$observedCounts, cellType = sim$z)
}

## Load a 10x PBMC data set and keep cells and genes with counts
loadPbmc <- function(dataset, file) {
  sce <- TENxPBMCData::TENxPBMCData(dataset = dataset)
  counts <- methods::as(SingleCellExperiment::counts(sce), "dgCMatrix")
  colnames(counts) <- sce$Barcode
  rownames(counts) <- SummarizedExperiment::rowData(sce)$ENSEMBL_ID
  counts <- counts[Matrix::rowSums(counts) > 0, Matrix::colSums(counts) > 0]
  saveRDS(list(counts = counts), file)
}

dataFile <- function(outdir, name) {
  dir.create(file.path(outdir, "data"), showWarnings = FALSE)
  file.path(outdir, "data", paste0(name, ".rds"))
}


## ---------------------------------------------------------------------
## Worker: runs one task in a fresh process
## ---------------------------------------------------------------------

## Current RSS of this process in MB (NA without the ps package)
currentRss <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) {
    return(NA_real_)
  }
  unname(ps::ps_memory_info()[["rss"]]) / 2^20
}

## Size of the CSV files rstan::vb() leaves in tempdir(), in MB
stanCsvMb <- function() {
  csv <- list.files(tempdir(), pattern = "\\.csv$", full.names = TRUE)
  sum(file.size(csv)) / 2^20
}

## EM iterations of each batch, as "n" or "n1;n2;..."
emIterations <- function(result) {
  if (is.null(result$estimates)) {
    return(NA_character_)
  }
  iters <- vapply(result$estimates, function(e) e$iteration, integer(1))
  paste(iters, collapse = ";")
}

## Clusters used, summed over batches
clusterCount <- function(result) {
  if (is.null(result$estimates)) {
    return(NA_integer_)
  }
  sum(vapply(result$estimates, function(e) nlevels(e$z), integer(1)))
}

## The call each decontX run type makes, given the loaded data
decontXCall <- function(run, data) {
  switch(run,
    z_supplied = function() {
      decontX::decontX(data$counts, z = data$z, verbose = FALSE)
    },
    z_estimated = function() {
      decontX::decontX(data$counts, verbose = FALSE)
    },
    batched4 = function() {
      decontX::decontX(data$counts,
        z = data$z, batch = data$batch,
        verbose = FALSE
      )
    },
    background = function() {
      decontX::decontX(data$counts,
        z = data$z, background = data$background,
        verbose = FALSE
      )
    },
    stop("Unknown run type: ", run)
  )
}

## Time `fn` with bench::mark and return one row of results. Memory is
## measured in a separate, untimed run (bench_memory() traces every R
## allocation, which slows the call), so a case runs iterations + 1
## times when `memory` is TRUE. cpu_s is the process CPU time per timed
## iteration; compare it with median_s to spot a throttled run.
markCase <- function(fn, iterations, memory = TRUE) {
  memAllocMb <- NA_real_
  if (memory) {
    traced <- bench::bench_memory(fn())
    memAllocMb <- as.numeric(traced$mem_alloc) / 2^20
    rm(traced)
    invisible(gc())
  }
  holder <- new.env()
  cpuStart <- proc.time()
  marked <- bench::mark(
    holder$result <- fn(),
    iterations = iterations,
    check = FALSE,
    memory = FALSE,
    filter_gc = FALSE
  )
  cpuUsed <- proc.time() - cpuStart
  times <- as.numeric(marked$time[[1]])
  list(
    row = data.frame(
      iterations = iterations,
      median_s = stats::median(times),
      min_s = min(times),
      max_s = max(times),
      cpu_s = sum(cpuUsed[c("user.self", "sys.self")]) / iterations,
      mem_alloc_mb = memAllocMb,
      n_gc = sum(marked$n_gc)
    ),
    result = holder$result
  )
}

## Run a decontX or decontPro timing case
workerBench <- function(task, outdir) {
  parts <- strsplit(task, "_", fixed = TRUE)[[1]]
  kind <- parts[1]
  label <- envOr("DECONTX_BENCH_LABEL", NA_character_)

  if (kind == "pro") {
    ## pro_<cells>
    cells <- as.integer(parts[2])
    adt <- simulateAdt(cells)
    rssBefore <- currentRss()
    csvBefore <- stanCsvMb()
    ## One run: Stan variational inference takes minutes, and its memory
    ## is allocated in C++, which bench_memory() cannot see
    marked <- markCase(function() {
      decontX::decontPro(adt$counts, adt$cellType)
    }, iterations = 1L, memory = FALSE)
    return(cbind(
      data.frame(
        label = label, mode = "sim", run = "decontPro", cells = cells,
        genes = nrow(adt$counts), nnz = sum(adt$counts > 0)
      ),
      marked$row,
      data.frame(
        rss_before_mb = rssBefore, em_iterations = NA_character_,
        clusters = proCellTypes, stan_csv_mb = stanCsvMb() - csvBefore
      )
    ))
  }

  ## sim_<cells>_<run> or real_<dataset>_<run>
  run <- paste(parts[-(1:2)], collapse = "_")
  if (kind == "sim") {
    cells <- as.integer(parts[2])
    data <- readRDS(dataFile(outdir, paste0("sim-", cells)))
  } else {
    data <- readRDS(dataFile(outdir, parts[2]))
    cells <- ncol(data$counts)
  }
  if (run == "sce_hdf5") {
    ## The 10x data as distributed: a SingleCellExperiment whose counts
    ## are an HDF5-backed DelayedMatrix (exercises class restoration)
    data$sce <- TENxPBMCData::TENxPBMCData(dataset = parts[2])
    fn <- function() decontX::decontX(data$sce, verbose = FALSE)
  } else {
    fn <- decontXCall(run, data)
  }
  invisible(gc())
  rssBefore <- currentRss()
  marked <- markCase(fn, iterationsFor(cells))
  result <- marked$result
  if (methods::is(result, "SingleCellExperiment")) {
    result <- S4Vectors::metadata(result)$decontX
  }
  cbind(
    data.frame(
      label = label, mode = kind, run = run, cells = cells,
      genes = nrow(data$counts), nnz = length(data$counts@x)
    ),
    marked$row,
    data.frame(
      rss_before_mb = rssBefore, em_iterations = emIterations(result),
      clusters = clusterCount(result), stan_csv_mb = NA_real_
    )
  )
}

## Profiling: classify each Rprof sample by the innermost known function
## on its call stack, so time splits into decontX phases. Compiled code
## is charged to the R function that calls it (for example decontXEM).
phaseRules <- c(
  fit_dirichlet = "em_fit_dirichlet",
  decontXEM = "em_decontXEM",
  decontXLogLik = "em_loglik",
  decontXInitialize = "em_initialize",
  calculateNativeMatrix = "native_matrix",
  run_umap = "init_umap_optimize",
  findKNN = "init_umap_knn",
  buildIndex = "init_umap_knn",
  runUmap = "init_umap_other",
  runPca = "init_pca",
  modelGeneVariances = "init_hvg",
  chooseHighlyVariableGenes = "init_hvg",
  normalizeCounts = "init_normalize",
  centerSizeFactors = "init_normalize",
  dbscan = "init_dbscan",
  kmeans = "init_kmeans",
  .decontxInitializeZ = "init_other",
  .decontXoneBatch = "batch_other",
  .decontX = "core_assembly_restore",
  `sampler$call_sampler` = "stan_sampler",
  read_one_stan_csv = "stan_read_csv",
  .call_stan_vb = "stan_rstan_post",
  .process_stan_vb_out = "stan_process_output",
  .decontPro = "decontPro_prep",
  decontX = "method_writeback",
  decontPro = "method_other"
)

## Read an Rprof file and return seconds per phase
profilePhases <- function(profFile) {
  lines <- readLines(profFile)
  header <- lines[1]
  interval <- as.numeric(sub(".*sample.interval=(\\d+).*", "\\1", header))
  samples <- lines[-1]
  samples <- samples[!grepl("^#", samples)]
  ## Drop the memory-profiling prefix (":n:n:n:n:") of each sample
  samples <- sub("^:[0-9:]*:", "", samples)
  stacks <- lapply(strsplit(samples, "\" \"", fixed = TRUE), function(s) {
    frames <- gsub("\"", "", trimws(s), fixed = TRUE)
    ## "pkg::fun" and "pkg:::fun" count as "fun"
    sub("^[A-Za-z0-9.]+:::?", "", frames)
  })
  phases <- vapply(stacks, function(stack) {
    known <- stack[stack %in% names(phaseRules)]
    if (length(known) == 0) "outside" else phaseRules[[known[1]]]
  }, character(1))
  seconds <- table(phases) * interval / 1e6
  out <- data.frame(
    phase = names(seconds),
    seconds = as.numeric(seconds)
  )
  out$percent <- round(100 * out$seconds / sum(out$seconds), 1)
  out[order(-out$seconds), ]
}

## Profile one call: Rprof, summary tables, phases, profvis HTML
profileCall <- function(name, fn, profDir) {
  profFile <- file.path(profDir, paste0(name, ".Rprof"))
  invisible(gc())
  rssBefore <- currentRss()
  utils::Rprof(profFile, interval = 0.01, memory.profiling = TRUE)
  elapsed <- system.time(result <- fn())[["elapsed"]]
  utils::Rprof(NULL)

  summary <- utils::summaryRprof(profFile, memory = "both")
  utils::write.csv(summary$by.total,
    file.path(profDir, paste0(name, "-by-total.csv"))
  )
  utils::write.csv(summary$by.self,
    file.path(profDir, paste0(name, "-by-self.csv"))
  )
  phases <- profilePhases(profFile)
  utils::write.csv(phases, file.path(profDir, paste0(name, "-phases.csv")),
    row.names = FALSE
  )

  ## profvis reads the same Rprof file; write HTML when it is available
  if (requireNamespace("profvis", quietly = TRUE)) {
    widget <- profvis::profvis(prof_input = profFile)
    htmlFile <- file.path(profDir, paste0(name, "-profvis.html"))
    try(htmlwidgets::saveWidget(widget, htmlFile, selfcontained = FALSE),
      silent = TRUE
    )
  }
  list(
    result = result, elapsed = elapsed, rss_before_mb = rssBefore,
    phases = phases
  )
}

## Run the profile task: profile_<cells>_<what>
workerProfile <- function(task, outdir) {
  parts <- strsplit(task, "_", fixed = TRUE)[[1]]
  cells <- as.integer(parts[2])
  what <- paste(parts[-(1:2)], collapse = "_")
  profDir <- file.path(outdir, "profile")
  dir.create(profDir, showWarnings = FALSE)

  if (what == "decontPro") {
    adt <- simulateAdt(cells)
    fn <- function() decontX::decontPro(adt$counts, adt$cellType)
  } else {
    ## A SingleCellExperiment, so SCE write-back is part of the profile
    data <- readRDS(dataFile(outdir, paste0("sim-", cells)))
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = data$counts)
    )
    z <- if (what == "z_estimated") NULL else data$z
    batch <- if (what == "batched4") data$batch else NULL
    fn <- function() {
      decontX::decontX(sce, z = z, batch = batch, verbose = FALSE)
    }
  }
  prof <- profileCall(task, fn, profDir)
  result <- prof$result
  if (methods::is(result, "SingleCellExperiment")) {
    result <- S4Vectors::metadata(result)$decontX
  }
  phases <- prof$phases
  data.frame(
    label = envOr("DECONTX_BENCH_LABEL", NA_character_),
    mode = "profile", run = what, cells = cells,
    elapsed_s = prof$elapsed, rss_before_mb = prof$rss_before_mb,
    em_iterations = emIterations(result),
    phases = paste0(phases$phase, "=", round(phases$seconds, 2),
      collapse = "; "
    )
  )
}

## Entry point of a worker process
runWorkerTask <- function(task, outdir, lib) {
  suppressPackageStartupMessages(
    library(decontX, lib.loc = lib)
  )
  loadedFrom <- dirname(find.package("decontX"))
  if (normalizePath(loadedFrom) != normalizePath(lib)) {
    stop("decontX loaded from ", loadedFrom, ", not ", lib)
  }
  if (!requireNamespace("bench", quietly = TRUE)) {
    stop("Install the 'bench' package to run the benchmark.")
  }

  if (startsWith(task, "simulate_")) {
    cells <- as.integer(sub("simulate_", "", task, fixed = TRUE))
    simulateData(cells, dataFile(outdir, paste0("sim-", cells)))
    result <- TRUE
  } else if (startsWith(task, "loadpbmc_")) {
    dataset <- sub("loadpbmc_", "", task, fixed = TRUE)
    loadPbmc(dataset, dataFile(outdir, dataset))
    result <- TRUE
  } else if (startsWith(task, "session_")) {
    sessionMode <- sub("session_", "", task, fixed = TRUE)
    writeLines(
      utils::capture.output(utils::sessionInfo()),
      file.path(outdir, paste0("sessionInfo-", sessionMode, ".txt"))
    )
    result <- TRUE
  } else if (startsWith(task, "profile_")) {
    result <- workerProfile(task, outdir)
  } else {
    result <- workerBench(task, outdir)
    result$decontx_version <- as.character(utils::packageVersion("decontX"))
  }
  saveRDS(result, file.path(outdir, "work", paste0(task, ".rds")))
}


## ---------------------------------------------------------------------
## Modes (parent process)
## ---------------------------------------------------------------------

## Results file for a mode, time-stamped so reruns never overwrite it
resultsFile <- function(outdir, mode) {
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  file.path(outdir, paste0("bench-", mode, "-", stamp, ".csv"))
}

## Run tasks one after another, binding the rows they return, and
## rewrite the CSV after each one so partial results survive
runTasks <- function(tasks, outdir, lib, csvFile) {
  rows <- list()
  for (task in tasks) {
    row <- runWorker(task, outdir, lib)
    if (is.data.frame(row)) {
      rows[[task]] <- row
      utils::write.csv(do.call(rbind, rows), csvFile, row.names = FALSE)
    }
  }
  if (length(rows) > 0) do.call(rbind, rows) else NULL
}

## Simulated data at each size, then decontPro
modeSim <- function(outdir, lib) {
  sizes <- as.integer(splitList(
    envOr("DECONTX_BENCH_SIZES", "10000,50000,100000")
  ))
  runs <- splitList(envOr(
    "DECONTX_BENCH_RUNS",
    "z_supplied,z_estimated,batched4,background"
  ))
  proCells <- as.integer(splitList(
    envOr("DECONTX_BENCH_PRO_CELLS", "2000,5000")
  ))
  csvFile <- resultsFile(outdir, "sim")

  tasks <- character()
  for (cells in sizes) {
    if (!file.exists(dataFile(outdir, paste0("sim-", cells)))) {
      runWorker(paste0("simulate_", cells), outdir, lib)
    }
    tasks <- c(tasks, paste("sim", cells, runs, sep = "_"))
  }
  if (length(proCells) > 0) {
    tasks <- c(tasks, paste("pro", proCells, sep = "_"))
  }
  results <- runTasks(tasks, outdir, lib, csvFile)
  runWorker("session_sim", outdir, lib)
  results
}

## Real PBMC data: dgCMatrix with estimated clusters, and the HDF5-backed
## SingleCellExperiment as distributed
modeReal <- function(outdir, lib) {
  dataset <- envOr("DECONTX_BENCH_PBMC", "pbmc4k")
  if (!file.exists(dataFile(outdir, dataset))) {
    runWorker(paste0("loadpbmc_", dataset), outdir, lib)
  }
  tasks <- paste("real", dataset, c("z_estimated", "sce_hdf5"), sep = "_")
  results <- runTasks(tasks, outdir, lib,
    resultsFile(outdir, "real")
  )
  runWorker("session_real", outdir, lib)
  results
}

## Rprof profiles of one mid-size decontX run (estimated clusters, user
## clusters, four batches with user clusters) and one decontPro run
modeProfile <- function(outdir, lib) {
  cells <- as.integer(envOr("DECONTX_BENCH_PROFILE_CELLS", "50000"))
  proCells <- as.integer(splitList(
    envOr("DECONTX_BENCH_PRO_CELLS", "2000,5000")
  ))
  if (!file.exists(dataFile(outdir, paste0("sim-", cells)))) {
    runWorker(paste0("simulate_", cells), outdir, lib)
  }
  tasks <- paste("profile", cells, c("z_estimated", "z_supplied", "batched4"),
    sep = "_"
  )
  if (length(proCells) > 0) {
    tasks <- c(tasks, paste("profile", max(proCells), "decontPro", sep = "_"))
  }
  results <- runTasks(tasks, outdir, lib,
    resultsFile(outdir, "profile")
  )
  runWorker("session_profile", outdir, lib)
  results
}


## ---------------------------------------------------------------------
## Main
## ---------------------------------------------------------------------

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  ## Worker processes: worker <task> <outdir> <lib>
  if (length(args) >= 1 && args[1] == "worker") {
    runWorkerTask(args[2], args[3], args[4])
    return(invisible())
  }

  mode <- if (length(args) >= 1) args[1] else "sim"
  if (!mode %in% c("sim", "real", "profile")) {
    stop("mode must be 'sim', 'real' or 'profile', not '", mode, "'.")
  }
  outdir <- if (length(args) >= 2) args[2] else tempfile("decontx-bench-")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outdir <- normalizePath(outdir)

  repoDir <- normalizePath(".")
  if (!file.exists(file.path(repoDir, "DESCRIPTION"))) {
    stop("Run this script from the decontX repository root.")
  }
  if (!nzchar(Sys.getenv("DECONTX_BENCH_LABEL"))) {
    commit <- tryCatch(
      system2("git",
        c("-C", shQuote(repoDir), "rev-parse", "--short", "HEAD"),
        stdout = TRUE, stderr = FALSE
      ),
      error = function(e) NA_character_
    )
    Sys.setenv(DECONTX_BENCH_LABEL = commit[1])
  }

  lib <- benchLibrary(repoDir, outdir)
  results <- switch(mode,
    sim = modeSim(outdir, lib),
    real = modeReal(outdir, lib),
    profile = modeProfile(outdir, lib)
  )
  message("Results in ", outdir)
  if (!is.null(results)) {
    print(results, row.names = FALSE)
  }
}

main()
