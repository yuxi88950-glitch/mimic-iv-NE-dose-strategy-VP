# Primary analysis: IPC-weighted KM, RMST and patient-cluster robust Cox.
run_analysis <- function(output_dir, B, workers, code_dir) {
  fitted <- fit_all(rep(1L, N), keep=TRUE, seed=20261017L)
  weights <- make_weights(fitted)
  absolute <- absolute_effects(weights$weight, keep=TRUE)
  cox <- coxph(Surv(tstart, tstop, event) ~ A, data=geom,
               weights=weights$weight, cluster=stay_id, robust=TRUE,
               ties='efron', control=coxph.control(timefix=FALSE))
  fwrite(absolute$curves, file.path(output_dir, 'survival.csv'))
  fwrite(fitted$monitor, file.path(output_dir, 'weight_models.csv'))
  fwrite(data.table(lower01=weights$lower, upper99=weights$upper),
         file.path(output_dir, 'weight_limits.csv'))

  # Frequency weights resample patients, keeping both clones together.
  # The numerator, ridge denominator, and truncation are refitted each time.
  cluster <- parallel::makePSOCKcluster(min(workers, B))
  on.exit(parallel::stopCluster(cluster), add=TRUE)
  parallel::clusterExport(cluster, 'code_dir', envir=environment())
  parallel::clusterCall(cluster, function(paths) .libPaths(paths), .libPaths())
  parallel::clusterEvalQ(cluster, {
    suppressPackageStartupMessages({
      library(data.table); library(Matrix); library(fastglm); library(glmnet)
    })
    setDTthreads(1)
    source(file.path(code_dir, 'R', 'weights.R'))
    NULL
  })
  parallel::clusterExport(cluster,
    c('N', 'panel', 'models', 'geom', 'clones', 'lastidx', 'gr', 'km_geometry'))
  bootstrap <- list()
  batches <- split(seq_len(B), ceiling(seq_len(B) / min(workers, B)))
  for (batch in batches) {
    answers <- parallel::parLapply(cluster, batch, function(b) {
      set.seed(20262000L + b)
      frequency <- tabulate(sample.int(N, N, replace=TRUE), nbins=N)
      fit <- fit_all(frequency, seed=20263000L + b)
      weight <- make_weights(fit, frequency)
      as.data.table(c(list(replicate=b),
        as.list(absolute_effects(weight$weight, frequency)$values)))
    })
    bootstrap <- c(bootstrap, answers)
    fwrite(rbindlist(bootstrap), file.path(output_dir, 'bootstrap.csv'))
    if (length(bootstrap) %% 20 == 0 || length(bootstrap) == B)
      cat(sprintf('Bootstrap: %d/%d\n', length(bootstrap), B))
    flush.console()
  }
  boot <- rbindlist(bootstrap)
  result <- rbindlist(lapply(c('risk_A25', 'risk_B35', 'RD',
                              'rmst_A25', 'rmst_B35', 'RMST_difference', 'HR'),
    function(estimand) {
      if (estimand == 'HR') {
        beta <- unname(coef(cox)[1]); se <- sqrt(vcov(cox)[1, 1])
        estimate <- exp(beta); ci <- exp(beta + c(-1, 1) * qnorm(.975) * se)
        method <- 'Patient-cluster robust Cox Wald'
      } else {
        estimate <- absolute$values[[estimand]]
        ci <- quantile(boot[[estimand]], c(.025, .975), names=FALSE)
        method <- sprintf('%d patient bootstrap percentiles', B)
      }
      data.table(estimand, estimate, lower95=ci[1], upper95=ci[2], method)
    }))
  fwrite(result, file.path(output_dir, 'primary.csv'))
  writeLines(c('Data: independently generated synthetic observations',
               paste('Patients:', N), paste('Clones:', nrow(clones)),
               paste('Bootstrap replicates:', B),
               'Contrasts: lower threshold minus/divided by higher threshold',
               'Risk estimates are proportions; RMST is measured in days.'),
             file.path(output_dir, 'run.txt'))
  draw_km(absolute$curves, output_dir)
  print(result)
  invisible(result)
}

draw_km <- function(curves, output_dir) {
  pdf(file.path(output_dir, 'KM.pdf'), width=7, height=5, useDingbats=FALSE)
  on.exit(dev.off())
  par(mar=c(4.3, 4.5, 2, 1))
  plot(0, 0, type='n', xlim=c(0, 28), ylim=c(0, 1),
       xlab='Days since time zero', ylab='IPC-weighted cumulative mortality',
       main='Synthetic-data demonstration', yaxt='n')
  axis(2, at=seq(0, 1, .2), labels=paste0(seq(0, 100, 20), '%'), las=1)
  for (arm_name in c('A25', 'B35')) {
    z <- curves[arm == arm_name]
    lines(z$time_h / 24, 1-z$survival, type='s', lwd=2,
          col=if (arm_name == 'A25') '#2166AC' else '#B35806')
  }
  legend('topleft', c('Lower threshold', 'Higher threshold'),
         col=c('#2166AC', '#B35806'), lwd=2, bty='n')
}
