# Rscript run.R [output_directory] [bootstrap=500] [workers=4]
if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', '.UTF-8'))
args <- commandArgs(TRUE)
file_arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)[1]
code_dir <- dirname(normalizePath(sub('^--file=', '', file_arg), winslash='/'))
output_dir <- if (length(args)) args[1] else file.path(code_dir, 'output')
B <- if (length(args) >= 2) as.integer(args[2]) else 500L
workers <- if (length(args) >= 3) as.integer(args[3]) else 4L
stopifnot(is.finite(B), B >= 1, is.finite(workers), workers >= 1)
packages <- c('data.table', 'survival', 'Matrix', 'fastglm', 'glmnet')
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop('Install required R packages: ', paste(missing, collapse=', '))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(Matrix)
  library(splines); library(fastglm); library(glmnet)
})
setDTthreads(1)
for (name in c('data.R', 'models.R', 'weights.R', 'analysis.R'))
  source(file.path(code_dir, 'R', name))
if (dir.exists(output_dir) && length(list.files(output_dir, all.files=TRUE, no..=TRUE)))
  stop('Choose a new or empty output directory.')
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
if (!dir.exists(output_dir)) stop('Cannot create output directory.')

cat('Synthetic-data demonstration; results are not study estimates.\n')
load_data(prepare_data(generate_data(n=800, seed=20261016)))
km_geometry <- build_km_geometry()
run_analysis(output_dir, B, workers, code_dir)
writeLines(capture.output(sessionInfo()), file.path(output_dir, 'sessionInfo.txt'))
cat('Completed:', normalizePath(output_dir, winslash='/'), '\n')
