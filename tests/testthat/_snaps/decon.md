# the bgBatch error and the cluster warning read correctly

    Code
      decontX(d$sim$observedCounts, z = d$sim$z, background = d$bg, bgBatch = d$
        bgBatch, verbose = FALSE)
    Condition
      Error in `.decontX()`:
      ! 'bgBatch' was given without 'batch'; give both or neither.

---

    Code
      sim <- simulateContamination(K = 101, C = 10)
    Condition
      Warning in `.simulateContaminatedMatrix()`:
      Only 10 clusters are simulated. Try to increase number of cells 'C' if more clusters are needed

