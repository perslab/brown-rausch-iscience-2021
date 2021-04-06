## Load your packages, e.g. library(drake).
source("./packages.R")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

## _drake.R must end with a call to drake_config().
## The arguments to drake_config() are basically the same as those to make().
## lock_envir allows functions that alter the random seed to be used. The biggest
## culprits of this seem to be interactive graphics e.g. plotly and mapdeck.
future::plan(future::multicore)
options(future.globals.maxSize= 50000*1024^2)
drake::drake_config(the_plan, parallelism = "future", 
                    jobs = 16, memory_strategy = "autoclean", 
                    format = "qs",
                    lock_envir = FALSE)

