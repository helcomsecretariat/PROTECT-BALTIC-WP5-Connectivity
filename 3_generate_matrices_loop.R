# Generate per-date matrices; run several sessions in parallel (run
# _matrix_run_check.R first). A failed date is logged and the loop continues.
for (j in 1:10) {
  msg <- tryCatch({ source("_generate_matrices.R"); "" },
                  error = function(e) conditionMessage(e))
  if (msg == "") next
  if (grepl("All matrices", msg)) { cat(msg, "\n"); break }   # nothing left to claim
  cat("Matrix run FAILED:", msg, "\n")
}
