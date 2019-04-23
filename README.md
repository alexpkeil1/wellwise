# wellwise: an R package

    devtools::install_github("alexpkeil1/wellwise")
    library(wellwise)
    ... 
    #stan/jags model code here
    ...
    # read all data
    RAW = data_importer(package='wellwise')
    # run the first  simulations
    res = analysis_wrapper(simiters=1, rawdata=RAW, dir="", root="example", 
          s.code=stanmodel, iter=4000, warmup=500, control=list(adapt_delta=0.8), 
          type='stan')

    ... 
    #do some diagnostics, e.g. (for stan):
    rr = As.mcmc.list(res[[1]])

    coda::effectiveSize(rr)
    BayesianTools::gelmanDiagnostics(BayesianTools::convertCoda(rr), plot=TRUE)
    ...
    # run the rest of the iterations after diagnostics look good
    res2 = analysis_wrapper(simiters=2:200, rawdata=RAW, dir="", root="horseshoe",
           s.code=stanmodel, stanfit = res[[1]])
    #summarize
    allres = do.call(c, list(res, res2))
    attributes(allres) <- attributes(res)
       



wellwise v0.5.21
