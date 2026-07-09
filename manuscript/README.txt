CCTP / SDDB manuscript package
==============================

Main deliverable
----------------
CCTP_SDBB_manuscript.pdf

Source
------
CCTP_SDBB_manuscript.md    Full English academic manuscript
references.bib             Bibliography
figures/                   Main and case-study figures
cctp_batch_summary.csv     N=10 results table (companion)

Author
------
Johel Padilla-Villanueva, DrPH
University of Puerto Rico, Medical Sciences Campus
ORCID: 0000-0002-5797-6931
https://github.com/johelpadilla/systemictau

Rebuild PDF (pandoc ≥ 2.x, xelatex, FreeSerif)
----------------------------------------------
/usr/local/opt/pandoc/bin/pandoc CCTP_SDBB_manuscript.md \
  -o CCTP_SDBB_manuscript.pdf \
  --pdf-engine=xelatex --citeproc \
  --resource-path=.:figures \
  -V geometry:margin=1in -V mainfont="FreeSerif" -V monofont="Menlo"
