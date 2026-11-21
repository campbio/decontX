# Decontaminate using decontPro

Decontaminate using decontPro

## Usage

``` r
decontPro(filtered_counts, cellType, ...)

# S4 method for class 'SingleCellExperiment'
decontPro(
  filtered_counts,
  cellType,
  deltaSd = 2e-05,
  backgroundSd = 2e-06,
  background = NULL,
  ...,
  iter = 50000,
  seed = 12345,
  verbose = TRUE
)

# S4 method for class 'Seurat'
decontPro(
  filtered_counts,
  cellType,
  deltaSd = 2e-05,
  backgroundSd = 2e-06,
  background = NULL,
  ...,
  iter = 50000,
  seed = 12345,
  verbose = TRUE
)

# S4 method for class 'ANY'
decontPro(
  filtered_counts,
  cellType,
  deltaSd = 2e-05,
  backgroundSd = 2e-06,
  background = NULL,
  ...,
  iter = 50000,
  seed = 12345,
  verbose = TRUE
)
```

## Arguments

- filtered_counts:

  Count matrix NxM (feature x droplet) with only filtered droplets after
  cell calling, with at least 2 ADTs and 2 droplets. The counts must be
  non-negative whole numbers with no NA. If this is a
  SingleCellExperiment or a Seurat object, expect `counts` in the assay
  slot. Split `counts` layers of a Seurat v5 assay are joined first, so
  the droplets follow `Cells()`.

- cellType:

  Vector of length M giving the cell type of each droplet, with no NA.
  Any labels can be used, such as integers, character strings or a
  factor. The rows of `r_est` in the returned parameters follow the
  sorted labels, or a factor's levels with unused levels dropped. Alias:
  `cell_type`.

- ...:

  Only the aliases `cell_type`, `delta_sd`, `background_sd` and
  `ambient_counts`, which can be used in place of `cellType`, `deltaSd`,
  `backgroundSd` and `background`. Giving both names of an argument, or
  any other argument, stops with an error.

- deltaSd:

  Standard deviation of the half-normal prior on how much the ambient
  contamination levels vary around each droplet's mean level. Smaller
  values keep the levels of a droplet closer together. Default 2e-5.
  Alias: `delta_sd`.

- backgroundSd:

  Standard deviation of the half-normal prior on how much the levels of
  the model's background component vary around each ADT's mean level.
  This background is part of the model, not the `background` counts.
  Default 2e-6. Alias: `background_sd`.

- background:

  Count matrix (feature x droplet) of empty droplets, used to estimate
  the ambient ADT profile; without it, the profile comes from
  `filtered_counts`. Only each ADT's total across these droplets (the
  row sums) is used, pooled into one ambient profile, so the droplets
  must be empty. Droplets whose barcodes (column names) are also in
  `filtered_counts` are removed, so the raw matrix with every barcode
  can be given. When either has no column names, nothing is removed,
  with a warning. A warning is also given when no barcode matches but
  `background` has at least a tenth as many droplets as there are cells
  with at least the median cell's total count, as a raw matrix does when
  the two name the same barcodes differently. Barcodes are matched
  across all droplets, so with several samples they must be unique
  across samples, for example with a sample prefix in both. Similar to
  `filtered_counts`, if it is a wrapper object, expect `counts` in the
  assay slot. When both have row names, its rows are matched to the ADTs
  of `filtered_counts` by name: ADTs that aren't in `filtered_counts`
  are dropped, and every ADT of `filtered_counts` must be present.
  Otherwise the rows are matched by position. The ADTs used must have no
  NA or negative counts, and some counts in total. Default NULL. Alias:
  `ambient_counts`.

- iter:

  Maximum number of iterations of Stan's variational inference. The fit
  usually stops earlier, when it converges; each iteration takes longer
  on larger data. Default 50000.

- seed:

  Seed for Stan's random number generator, a whole number from 0 to
  `.Machine$integer.max`. Default 12345.

- verbose:

  Whether to print messages and Stan's progress. With `FALSE`, Stan's
  output is shown only if the fit fails. Default TRUE.

## Value

A list. Each count of `filtered_counts` is split into three parts, which
add up to it:

- decontaminated_counts:

  The native part: the counts from the cell's own proteins. An ADT x
  droplet matrix with the row and column names of `filtered_counts`.

- ambient_counts:

  The ambient part, a matrix like `decontaminated_counts`.

- background_counts:

  The part from the model's background component, a matrix like
  `decontaminated_counts`.

- parameters:

  A list of the model's estimates: `delta_sd` and `background_sd`, the
  values of `deltaSd` and `backgroundSd` used; `r_est`, the native ADT
  profile of each cell type (cell type x ADT, rows named by the
  `cellType` labels, each row summing to 1); `p_est`, the ambient ADT
  profile (one value per ADT, summing to 1, named by the row names of
  the matrix it comes from, if any); `delta_est`, the weight of the
  ambient profile in the non-background part of each ADT in each
  droplet, the native profile getting the rest (ADT x droplet);
  `delta_mean_est`, its mean level in each droplet; `background_est`,
  the background as a proportion of the droplet's total count (ADT x
  droplet); and `background_mean_est`, its mean level for each ADT. The
  ADTs and droplets of these estimates are named by the row and column
  names of `filtered_counts`, if it has them.

## Examples

``` r
# Simulated counts of 10 ADTs (rows) in 8 droplets (columns)
set.seed(1)
counts <- matrix(rpois(80, lambda = 20), nrow = 10)

# Name the droplets by their barcodes
colnames(counts) <- paste0("cell", 1:8)

# The cell type of each droplet
cellTypes <- c("A", "A", "A", "A", "B", "B", "B", "B")

# Simulated counts of the same 10 ADTs in 8 empty droplets (optional)
emptyDroplets <- matrix(rpois(80, lambda = 2), nrow = 10)
colnames(emptyDroplets) <- paste0("empty", 1:8)

# Decontaminate the counts
out <- decontPro(counts, cellTypes, background = emptyDroplets)
#> Chain 1: ------------------------------------------------------------
#> Chain 1: EXPERIMENTAL ALGORITHM:
#> Chain 1:   This procedure has not been thoroughly tested and may be unstable
#> Chain 1:   or buggy. The interface is subject to change.
#> Chain 1: ------------------------------------------------------------
#> Chain 1: 
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Gradient evaluation took 7.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.72 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Begin eta adaptation.
#> Chain 1: Iteration:   1 / 250 [  0%]  (Adaptation)
#> Chain 1: Iteration:  50 / 250 [ 20%]  (Adaptation)
#> Chain 1: Iteration: 100 / 250 [ 40%]  (Adaptation)
#> Chain 1: Iteration: 150 / 250 [ 60%]  (Adaptation)
#> Chain 1: Iteration: 200 / 250 [ 80%]  (Adaptation)
#> Chain 1: Success! Found best value [eta = 1] earlier than expected.
#> Chain 1: 
#> Chain 1: Begin stochastic gradient ascent.
#> Chain 1:   iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes 
#> Chain 1:    100  -4992173899.365             1.000            1.000
#> Chain 1:    200    -63748884.200            39.155           77.310
#> Chain 1:    300      -902095.267            49.326           69.668
#> Chain 1:    400       -13437.765            53.527           69.668
#> Chain 1:    500        -1947.551            44.002           66.131
#> Chain 1:    600        -1601.295            36.704           66.131
#> Chain 1:    700        -1699.128            31.469            5.900
#> Chain 1:    800        -1120.261            27.600            5.900
#> Chain 1:    900        -1370.414            24.554            1.000
#> Chain 1:   1000        -1245.820            22.108            1.000
#> Chain 1:   1100        -1168.730            20.104            0.517   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1200        -1089.606            18.435            0.517   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1300        -1031.840            17.021            0.216   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1400        -1310.076            15.821            0.216   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1500         -902.802            14.796            0.216   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1600         -913.231            13.872            0.216   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1700         -976.686            13.060            0.212   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1800         -951.063            12.336            0.212   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1900         -815.643            11.695            0.183   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2000         -896.875            11.115            0.183   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2100         -946.521            10.588            0.166   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2200         -780.049            10.117            0.183   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2300         -723.013             9.680            0.166   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2400         -880.988             9.284            0.179   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2500         -737.759             8.921            0.179   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2600         -818.701             8.581            0.179   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2700         -776.930             8.266            0.166   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2800         -740.912             7.972            0.166   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2900         -779.084             7.699            0.100   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3000         -803.316             7.443            0.100   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3100         -793.438             7.204            0.099   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3200         -773.454             6.979            0.099   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3300         -810.002             6.769            0.091   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3400         -766.170             6.572            0.091   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3500         -683.048             6.387            0.091   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3600         -729.210             6.212            0.091   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3700         -727.476             6.044            0.079   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3800         -772.236             5.886            0.079   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3900         -719.570             5.737            0.073   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4000         -687.423             5.595            0.073   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4100         -690.571             5.459            0.073   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4200         -713.881             5.330            0.073   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4300         -679.750             5.207            0.066   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4400         -684.765             5.089            0.066   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4500         -646.449             4.977            0.065   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4600         -698.921             4.870            0.066   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4700         -663.483             4.768            0.065   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4800         -702.472             4.670            0.065   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4900         -669.846             4.575            0.063   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5000         -728.250             4.485            0.065   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5100         -656.841             4.468            0.065   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5200         -684.848             2.922            0.063   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5300         -692.807             1.529            0.059   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5400         -702.401             0.207            0.058
#> Chain 1:   5500         -676.205             0.090            0.058
#> Chain 1:   5600         -642.302             0.086            0.057
#> Chain 1:   5700         -620.167             0.086            0.056
#> Chain 1:   5800         -693.809             0.078            0.056
#> Chain 1:   5900         -680.708             0.074            0.056
#> Chain 1:   6000         -637.959             0.074            0.056
#> Chain 1:   6100         -659.858             0.073            0.054
#> Chain 1:   6200         -687.317             0.072            0.053
#> Chain 1:   6300         -680.027             0.071            0.053
#> Chain 1:   6400         -673.267             0.067            0.052
#> Chain 1:   6500         -642.839             0.059            0.050
#> Chain 1:   6600         -657.524             0.060            0.050
#> Chain 1:   6700         -651.914             0.058            0.049
#> Chain 1:   6800         -614.753             0.059            0.050
#> Chain 1:   6900         -647.156             0.057            0.050
#> Chain 1:   7000         -619.487             0.056            0.049
#> Chain 1:   7100         -618.079             0.055            0.049
#> Chain 1:   7200         -612.702             0.051            0.049
#> Chain 1:   7300         -635.092             0.050            0.047
#> Chain 1:   7400         -637.807             0.046            0.047
#> Chain 1:   7500         -613.619             0.043            0.045
#> Chain 1:   7600         -605.355             0.042            0.045
#> Chain 1:   7700         -651.336             0.042            0.045
#> Chain 1:   7800         -626.497             0.042            0.041
#> Chain 1:   7900         -606.074             0.041            0.040
#> Chain 1:   8000         -609.562             0.041            0.040
#> Chain 1:   8100         -626.087             0.041            0.040
#> Chain 1:   8200         -622.721             0.041            0.040
#> Chain 1:   8300         -609.909             0.040            0.040
#> Chain 1:   8400         -610.870             0.039            0.039
#> Chain 1:   8500         -609.821             0.037            0.039
#> Chain 1:   8600         -666.629             0.037            0.039
#> Chain 1:   8700         -621.851             0.039            0.039
#> Chain 1:   8800         -604.195             0.038            0.039
#> Chain 1:   8900         -619.582             0.037            0.036
#> Chain 1:   9000         -603.877             0.037            0.035
#> Chain 1:   9100         -596.097             0.037            0.035
#> Chain 1:   9200         -609.540             0.037            0.035
#> Chain 1:   9300         -613.718             0.036            0.034
#> Chain 1:   9400         -586.617             0.037            0.035
#> Chain 1:   9500         -586.953             0.035            0.034
#> Chain 1:   9600         -603.210             0.034            0.033
#> Chain 1:   9700         -603.220             0.033            0.029
#> Chain 1:   9800         -579.695             0.033            0.029
#> Chain 1:   9900         -609.793             0.033            0.029
#> Chain 1:   10000         -575.462             0.033            0.029
#> Chain 1:   10100         -614.677             0.032            0.029
#> Chain 1:   10200         -600.856             0.031            0.027
#> Chain 1:   10300         -585.681             0.032            0.027
#> Chain 1:   10400         -600.116             0.032            0.027
#> Chain 1:   10500         -593.953             0.031            0.026
#> Chain 1:   10600         -600.260             0.031            0.026
#> Chain 1:   10700         -571.321             0.031            0.026
#> Chain 1:   10800         -623.317             0.030            0.026
#> Chain 1:   10900         -599.025             0.031            0.026
#> Chain 1:   11000         -597.816             0.029            0.026
#> Chain 1:   11100         -618.144             0.029            0.026
#> Chain 1:   11200         -586.523             0.030            0.026
#> Chain 1:   11300         -569.370             0.030            0.026
#> Chain 1:   11400         -588.178             0.031            0.027
#> Chain 1:   11500         -567.784             0.030            0.027
#> Chain 1:   11600         -602.370             0.031            0.029
#> Chain 1:   11700         -575.140             0.032            0.030
#> Chain 1:   11800         -608.117             0.032            0.030
#> Chain 1:   11900         -569.670             0.032            0.030
#> Chain 1:   12000         -607.955             0.032            0.030
#> Chain 1:   12100         -592.390             0.033            0.030
#> Chain 1:   12200         -577.624             0.033            0.030
#> Chain 1:   12300         -600.525             0.033            0.030
#> Chain 1:   12400         -590.159             0.034            0.030
#> Chain 1:   12500         -574.558             0.033            0.029
#> Chain 1:   12600         -599.026             0.034            0.030
#> Chain 1:   12700         -603.806             0.033            0.029
#> Chain 1:   12800         -583.032             0.033            0.029
#> Chain 1:   12900         -579.357             0.032            0.027
#> Chain 1:   13000         -568.408             0.032            0.027
#> Chain 1:   13100         -568.704             0.032            0.027
#> Chain 1:   13200         -581.177             0.032            0.027
#> Chain 1:   13300         -576.542             0.032            0.027
#> Chain 1:   13400         -566.903             0.032            0.027
#> Chain 1:   13500         -583.232             0.033            0.028
#> Chain 1:   13600         -581.999             0.031            0.027
#> Chain 1:   13700         -575.379             0.030            0.027
#> Chain 1:   13800         -616.791             0.031            0.027
#> Chain 1:   13900         -596.219             0.031            0.027
#> Chain 1:   14000         -596.892             0.030            0.027
#> Chain 1:   14100         -576.327             0.031            0.028
#> Chain 1:   14200         -605.126             0.031            0.030
#> Chain 1:   14300         -576.454             0.032            0.032
#> Chain 1:   14400         -577.162             0.031            0.030
#> Chain 1:   14500         -577.448             0.031            0.030
#> Chain 1:   14600         -578.727             0.031            0.030
#> Chain 1:   14700         -564.857             0.031            0.030
#> Chain 1:   14800         -580.122             0.031            0.028
#> Chain 1:   14900         -570.423             0.030            0.027
#> Chain 1:   15000         -579.022             0.029            0.026
#> Chain 1:   15100         -568.635             0.028            0.026
#> Chain 1:   15200         -565.732             0.028            0.026
#> Chain 1:   15300         -571.210             0.028            0.026
#> Chain 1:   15400         -575.916             0.027            0.026
#> Chain 1:   15500         -559.591             0.028            0.026
#> Chain 1:   15600         -570.692             0.028            0.026
#> Chain 1:   15700         -565.756             0.027            0.026
#> Chain 1:   15800         -608.992             0.027            0.026
#> Chain 1:   15900         -559.845             0.028            0.026
#> Chain 1:   16000         -573.594             0.028            0.026
#> Chain 1:   16100         -559.529             0.028            0.026
#> Chain 1:   16200         -586.645             0.028            0.026
#> Chain 1:   16300         -578.710             0.028            0.025
#> Chain 1:   16400         -573.014             0.027            0.025
#> Chain 1:   16500         -572.726             0.026            0.024
#> Chain 1:   16600         -569.479             0.025            0.021
#> Chain 1:   16700         -577.361             0.025            0.019
#> Chain 1:   16800         -611.945             0.025            0.019
#> Chain 1:   16900         -563.848             0.025            0.019
#> Chain 1:   17000         -581.347             0.024            0.019
#> Chain 1:   17100         -584.380             0.024            0.019
#> Chain 1:   17200         -568.823             0.024            0.019
#> Chain 1:   17300         -563.622             0.024            0.018
#> Chain 1:   17400         -557.267             0.023            0.018
#> Chain 1:   17500         -570.753             0.023            0.018
#> Chain 1:   17600         -557.103             0.023            0.018
#> Chain 1:   17700         -581.970             0.024            0.019
#> Chain 1:   17800         -559.338             0.024            0.019
#> Chain 1:   17900         -560.803             0.024            0.019
#> Chain 1:   18000         -574.321             0.024            0.019
#> Chain 1:   18100         -556.779             0.024            0.021
#> Chain 1:   18200         -554.526             0.024            0.019
#> Chain 1:   18300         -556.651             0.024            0.019
#> Chain 1:   18400         -559.861             0.024            0.019
#> Chain 1:   18500         -592.543             0.024            0.019
#> Chain 1:   18600         -557.815             0.026            0.024
#> Chain 1:   18700         -570.362             0.026            0.024
#> Chain 1:   18800         -566.537             0.025            0.022
#> Chain 1:   18900         -576.224             0.024            0.019
#> Chain 1:   19000         -561.421             0.025            0.022
#> Chain 1:   19100         -561.534             0.024            0.019
#> Chain 1:   19200         -563.441             0.023            0.018
#> Chain 1:   19300         -569.466             0.022            0.017
#> Chain 1:   19400         -572.938             0.022            0.017
#> Chain 1:   19500         -552.301             0.023            0.018
#> Chain 1:   19600         -572.738             0.024            0.019
#> Chain 1:   19700         -556.949             0.024            0.019
#> Chain 1:   19800         -551.520             0.024            0.018
#> Chain 1:   19900         -554.380             0.023            0.018
#> Chain 1:   20000         -560.019             0.023            0.018
#> Chain 1:   20100         -559.309             0.023            0.017
#> Chain 1:   20200         -562.261             0.023            0.017
#> Chain 1:   20300         -558.810             0.023            0.017
#> Chain 1:   20400         -585.759             0.024            0.019
#> Chain 1:   20500         -552.283             0.024            0.019
#> Chain 1:   20600         -550.948             0.024            0.017
#> Chain 1:   20700         -553.740             0.024            0.017
#> Chain 1:   20800         -575.396             0.023            0.017
#> Chain 1:   20900         -563.558             0.022            0.017
#> Chain 1:   21000         -564.493             0.021            0.014
#> Chain 1:   21100         -559.835             0.021            0.014
#> Chain 1:   21200         -559.735             0.020            0.011
#> Chain 1:   21300         -562.791             0.020            0.011
#> Chain 1:   21400         -561.415             0.020            0.011
#> Chain 1:   21500         -551.692             0.020            0.011
#> Chain 1:   21600         -552.700             0.020            0.011
#> Chain 1:   21700         -562.982             0.020            0.011
#> Chain 1:   21800         -561.060             0.019            0.011
#> Chain 1:   21900         -567.537             0.018            0.011
#> Chain 1:   22000         -554.569             0.017            0.011
#> Chain 1:   22100         -556.495             0.017            0.011
#> Chain 1:   22200         -571.697             0.017            0.011
#> Chain 1:   22300         -551.669             0.018            0.011
#> Chain 1:   22400         -566.824             0.018            0.011
#> Chain 1:   22500         -554.952             0.018            0.011
#> Chain 1:   22600         -547.594             0.018            0.011
#> Chain 1:   22700         -554.126             0.017            0.011
#> Chain 1:   22800         -551.920             0.017            0.011
#> Chain 1:   22900         -547.394             0.017            0.011
#> Chain 1:   23000         -562.239             0.017            0.011
#> Chain 1:   23100         -551.108             0.017            0.011
#> Chain 1:   23200         -543.051             0.017            0.011
#> Chain 1:   23300         -562.501             0.017            0.012
#> Chain 1:   23400         -555.106             0.018            0.013
#> Chain 1:   23500         -544.922             0.017            0.013
#> Chain 1:   23600         -561.692             0.016            0.013
#> Chain 1:   23700         -551.535             0.016            0.013
#> Chain 1:   23800         -551.951             0.016            0.013
#> Chain 1:   23900         -555.155             0.016            0.012
#> Chain 1:   24000         -551.379             0.015            0.011
#> Chain 1:   24100         -555.083             0.015            0.011
#> Chain 1:   24200         -552.014             0.016            0.011
#> Chain 1:   24300         -542.835             0.016            0.012
#> Chain 1:   24400         -546.505             0.016            0.012
#> Chain 1:   24500         -557.004             0.015            0.012
#> Chain 1:   24600         -558.942             0.015            0.011
#> Chain 1:   24700         -555.873             0.014            0.010
#> Chain 1:   24800         -552.121             0.014            0.010
#> Chain 1:   24900         -550.136             0.014            0.010
#> Chain 1:   25000         -569.522             0.015            0.011
#> Chain 1:   25100         -565.896             0.015            0.011
#> Chain 1:   25200         -564.454             0.015            0.011
#> Chain 1:   25300         -566.075             0.015            0.011
#> Chain 1:   25400         -553.137             0.014            0.011
#> Chain 1:   25500         -548.490             0.013            0.008   MEDIAN ELBO CONVERGED
#> Chain 1: 
#> Chain 1: Drawing a sample of size 2 from the approximate posterior... 
#> Chain 1: COMPLETED.

# The decontaminated counts
decontaminatedCounts <- out$decontaminated_counts
```
