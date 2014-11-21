# Analysis of sieve effects in the RTS,S malaria vaccine trial

## Data prep

The haplotype data files from the Broad are:

```
qdata/sequences/BEP.tsv
qdata/sequences/CST3.tsv
qdata/sequences/SERA2.tsv
qdata/sequences/TEP.tsv
qdata/sequences/Th2R_Th3R.tsv
qdata/sequences/Th2R.tsv
qdata/sequences/Th3R.tsv
qdata/sequences/TRAP.tsv
qdata/sequences/Unnamed.tsv
```

Necessary clinical data is at:

```
qdata/clinical/sample_data.tsv
adata/RTSSclinicalData.csv
```

Together, these are compiled to:

```
adata/marks_data_c.tsv
adata/marks_data_x.tsv
```

with

```
python scripts/prep-genotype-data.py C > adata/marks_data_c.tsv
python scripts/prep-genotype-data.py X > adata/marks_data_x.tsv
```

## Descriptive statistics of genotype data

* [Results of genotype analysis](descriptive-analysis/descriptive-analysis.md)
