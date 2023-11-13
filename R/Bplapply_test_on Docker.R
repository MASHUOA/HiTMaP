fun <- function(v) {
  sqrt(v)
}

BPPARAM2=HiTMaP:::Parallel.OS(2)
BPPARAM4=HiTMaP:::Parallel.OS(4)
BPPARAM8=HiTMaP:::Parallel.OS(8)

start.time <- Sys.time()
lapply(1:1000000, fun) -> res
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken -> res

start.time <- Sys.time()
bplapply(1:1000000, fun,BPPARAM = BPPARAM2) -> res2
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken -> res2

BPPARAM2f=BPPARAM2
BPPARAM2f$fallback=F
BPPARAM2f$stop.on.error=F

start.time <- Sys.time()
bplapply(1:1000000, fun,BPPARAM = BPPARAM2f) -> res2f
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken -> res2f

start.time <- Sys.time()
bplapply(1:1000000, fun,BPPARAM = BPPARAM4) -> res4
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken -> res4


start.time <- Sys.time()
start.time <- Sys.time()
bplapply(1:1000000, fun,BPPARAM = BPPARAM8) -> res8
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken -> res8




