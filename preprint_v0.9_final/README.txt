CCTP v0.9_final Preprint Package
=================================

This directory contains the preprint-ready materials (v0.9) for:

  Context-Dependent Relational Reorganization of Heart Rate Dynamics
  Precedes Spontaneous Ventricular Fibrillation: Systemic Tau and Ordinal
  RECD Evidence from the Sudden Cardiac Death Holter Database

Contents:
- CCTP_SDBB_preprint_draft.md   (main manuscript draft)
- cctp_batch_summary.csv        (N=10 consolidated results)
- figures/batch/                (comparative batch plots)
- figures/{30,35,38,50}/        (key per-record diagnostic figures, including fixed 06_ews_panels.png)

Re-calibration used:
  theta3=0.08, high_thresh=0.65, lambda-relative

PDF (approved v0.9 ready):
  CCTP_SDBB_preprint_v0.9.pdf   (5 pages, FreeSerif / xelatex, Unicode-safe)

To re-render the draft locally (use pandoc ≥ 2.x / 3.x, not the macOS 1.13 default):
  /usr/local/opt/pandoc/bin/pandoc CCTP_SDBB_preprint_draft.md \
    -o CCTP_SDBB_preprint_v0.9.pdf --pdf-engine=xelatex \
    -V geometry:margin=1in -V mainfont="FreeSerif" -V monofont="Menlo"

Reproducibility commands (final parameters):
  python3 code/run_cctp_batch.py \
    --records 30,31,32,35,36,38,45,47,50,51 \
    --theta3 0.08 --high-thresh 0.65 --lambda-relative --force

Status: Preprint draft v0.9 + re-calibration completed. Ready for medRxiv / Preprints.org upload.

Generated: 2026-07-08
