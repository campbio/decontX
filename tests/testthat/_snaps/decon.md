# misspelled messages keep their current wording

    Code
      decontX(d$sim$observedCounts, z = d$sim$z, background = d$bg, bgBatch = d$
        bgBatch, verbose = FALSE)
    Condition
      Error in `.decontX()`:
      ! When experiment default to no bacth, background should also default to no batch.

---

    Code
      sim <- simulateContamination(K = 101, C = 10)
    Condition
      Warning in `.simulateContaminatedMatrix()`:
      Only 10 clusters are simulated. Try to increase numebr of cells 'C' if more clusters are needed

