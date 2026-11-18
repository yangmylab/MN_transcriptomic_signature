# Environment record

Generated: 2026-10-19 18:26:30 CST
R: R version 4.6.1 (2026-06-24) on x86_64-pc-linux-gnu

## Rebuild on another machine

1. Install R 4.6.x and Bioconductor 3.23.
2. Restore the repository packages from the lockfile:

   ```r
   install.packages("renv")
   renv::restore(lockfile = "environment/renv.lock")
   ```

3. The locally built Brainarray CDF packages are NOT in any repository. Rebuild them from the
   `.cdf.gz` files inside the GEO archives (no network needed):

   ```
   Rscript R/00d_install_cdf.R      # makecdfenv into ./.Rlib
   Rscript R/02a_build_cdf.R        # build the 5 CDF source packages from $DATA/*/raw/
   for p in work/cdf/*; do R CMD INSTALL -l .Rlib "$p"; done
   ```

4. Point the data root at the shared dataset collection (or copy it next to the project):

   ```
   export PAPER1_DATA=/path/to/MN/Datasets
   ```

5. Verify: all six verifiers must exit 0.

   ```
   for v in phase_a phase_b downstream scrna phase_d phase_e; do Rscript R/99_verify_$v.R; done
   ```

## Notes

- The analysis library is `./.Rlib` (prepended by `R/_lib.R`); the user library is
  `~/R/x86_64-pc-linux-gnu-library/4.6`. Both are recorded in `packages_*.csv`.
- `renv` is intentionally NOT activated at run time; `renv.lock` is a reproducibility record.
- Packages with `Source = Local` are the Brainarray CDF builds listed above.
