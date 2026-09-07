# Original data — not included

The files analysed here were released with the original Article and are **not
redistributed in this repository**. Download them from:

**Code Ocean capsule 9843791**
(linked from Ashokkumar, Hewitt, Ghezae & Willer, *Nature*, 2026,
https://doi.org/10.1038/s41586-026-10742-x)

Place the following in this directory:

- `rct_responses.RDS`
- `llm_responses.RDS`
- `forecasting_responses.RDS`
- `megastudies.RDS`

Then run `Rscript run_all.R` from the project root.

If you only want to reproduce the reported numbers, you do not need these
files: the derived comparison-level datasets in `data/derived/` are deposited
with this repository, and `Rscript run_all.R --from-derived` reproduces every
quantity in the manuscript from them alone.
