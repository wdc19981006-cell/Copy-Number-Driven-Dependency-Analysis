# R environment validation

R 4.5.0 is required. Renv is project-local; cache, sandbox and library paths are inside the project. Windows UTF-8 LC_CTYPE supports the en dash in the project name. Required package versions are in data/manifests/R_package_versions.csv and renv.lock.

The execution environment initially omitted PROCESSOR_ARCHITECTURE. Native cli shutdown called strcmp(getenv("PROCESSOR_ARCHITECTURE"), "ARM64") with NULL, causing Windows access violation 0xc0000005 in ucrtbase.dll after valid R outputs. Minimal probes and an owned-process debugger isolated the fault. Recompilation and a clean official R runtime did not remove the environment cause. Windows GetNativeSystemInfo returned architecture code 9 (AMD64). .Rprofile and the analysis entry point now fill this missing variable from the verified executable architecture before packages load. Original locked cli/rlang versions were restored; the final probe exited normally.

Some current CRAN Windows binaries report that they were built under R 4.5.3; the actual interpreter remains 4.5.0. Successful package loading, actual analysis exit codes and renv::status provide runtime validation. Earlier failed native probes/rebuilds are diagnostic history, not successful case runs.
