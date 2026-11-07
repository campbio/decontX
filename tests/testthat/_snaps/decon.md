# decontX writes the same log text

    Code
      writeLines(maskLog(readLines(logfile)))
    Output
      -------------------------------------------------- |
      Starting DecontX |
      -------------------------------------------------- |
      <date> .. Analyzing all cells |
      <date> .... Converting to sparse matrix |
      <date> .... Generating UMAP |
      <date> .... Estimating contamination |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> .. Calculating final decontaminated matrix |
      -------------------------------------------------- |
      Completed DecontX. Total time: <time> |
      -------------------------------------------------- |

---

    Code
      writeLines(maskLog(readLines(logfile)))
    Output
      <date> ..  1  cells in the background matrix were removed as they were found in  the filtered matrix. |
      -------------------------------------------------- |
      Starting DecontX |
      -------------------------------------------------- |
      <date> .. Dropped genes in 'background' that aren't in 'x': 1 of 61 (extra)|
      <date> .. Analyzing cells in batch 'A'|
      <date> .... Converting to sparse matrix |
      <date> .... Generating UMAP |
      <date> .... Estimating contamination |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> .. Calculating final decontaminated matrix |
      <date> .. Analyzing cells in batch 'B'|
      <date> .... Converting to sparse matrix |
      <date> .... Generating UMAP |
      <date> .... Estimating contamination |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> .. Calculating final decontaminated matrix |
      <date> .. Converting decontaminated matrix to a DelayedMatrix |
      -------------------------------------------------- |
      Completed DecontX. Total time: <time> |
      -------------------------------------------------- |

---

    Code
      writeLines(maskLog(readLines(logfile)))
    Output
      -------------------------------------------------- |
      Starting DecontX |
      -------------------------------------------------- |
      <date> .. Analyzing all cells |
      <date> .... Converting to sparse matrix |
      <date> .... Generating UMAP and estimating cell types |
      <date> .... Estimating contamination |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> .. Calculating final decontaminated matrix |
      -------------------------------------------------- |
      Completed DecontX. Total time: <time> |
      -------------------------------------------------- |

# decontX writes the same messages to the console

    Code
      writeLines(maskLog(sub("\n$", "", messages)))
    Output
      --------------------------------------------------|
      Starting DecontX|
      --------------------------------------------------|
      <date> .. Analyzing all cells|
      <date> .... Converting to sparse matrix|
      <date> .... Generating UMAP|
      <date> .... Estimating contamination|
      <date> ...... Completed iteration: <n> | converge: <x>|
      <date> ...... Completed iteration: <n> | converge: <x>|
      <date> ...... Completed iteration: <n> | converge: <x>|
      <date> .. Calculating final decontaminated matrix|
      --------------------------------------------------|
      Completed DecontX. Total time: <time>|
      --------------------------------------------------|

# the bgBatch error and the cluster warning read correctly

    Code
      decontX(d$sim$observedCounts, z = d$sim$z, background = d$bg, bgBatch = d$
        bgBatch, verbose = FALSE)
    Condition
      Error in `.decontXBatches()`:
      ! 'bgBatch' was given without 'batch'; give both or neither.

---

    Code
      sim <- simulateContamination(K = 101, C = 10)
    Condition
      Warning in `.simulateContaminatedMatrix()`:
      Only 10 clusters are simulated. Try to increase number of cells 'C' if more clusters are needed

# embedding = FALSE leaves the UMAP out of the log

    Code
      writeLines(maskLog(readLines(logfile)))
    Output
      -------------------------------------------------- |
      Starting DecontX |
      -------------------------------------------------- |
      <date> .. Analyzing all cells |
      <date> .... Converting to sparse matrix |
      <date> .... Estimating contamination |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> ...... Completed iteration: <n> | converge: <x> |
      <date> .. Calculating final decontaminated matrix |
      -------------------------------------------------- |
      Completed DecontX. Total time: <time> |
      -------------------------------------------------- |

