# Analysis of sieve effects in the RTS,S malaria vaccine trial

## Endpoints

1. Clinical malaria at 14 months
2. Severe clinical malaria at 14 months
3. Parasite infection at 20 months

## Loci

### CSP regions

1. TEP
2. CST3
3. Th2R (17 residues in length, starts at site 311)
4. Th3R (12 residues in length, starts at site 352)
5. Unnamed (12 residues in length, starts at site 294)
6. Th2R+Th3R Linkage

### CSP repeats

7. BEP

### Control regions

8. SERA2
9. TRAP

## Study sites

### Major sites

1. Agogo, Gabon, West Africa (41092)
2. Kintampo, Ghana, West Africa (41089)
3. Kombewa, Kenya, East Africa (41108)
4. Nanoro, Burkina Faso, West Africa (53530)
5. Siaya, Kenya, East Africa (53598)

### Minor sites

6. Bagamoyo, Tanzania, East Africa (41106)
7. Kilifi, Kenya, East Africa (41104)
8. Korogwe, Tanzania, East Africa (41107)
9. Lambarene, Gabon, West Africa (41087)
10. Lilongwe, Malawi, East Africa (53597)
11. Manhica, Mozambique, East Africa (41105)

## Data prep

### Marks data

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

### Site-specific marks data

Amended sequence data files are generated with:

```
python scripts/append-site-specific-match.py TEP > qdata/sequences/TEP_sites.tsv
python scripts/append-site-specific-match.py SERA2 > qdata/sequences/SERA2_sites.tsv
```

And site-specific mark files are then generated with:

```
python scripts/prep-genotype-data-sites.py C > adata/marks_data_c_sites.tsv
python scripts/prep-genotype-data-sites.py X > adata/marks_data_x_sites.tsv
```


## Descriptive statistics of genotype data

* [Code for genotype analysis](descriptive-analysis/descriptive-analysis.Rmd)

Compile to Markdown in R with:

```r
library(knitr)
knit("descriptive-analysis.Rmd")
```

* [Results of genotype analysis](descriptive-analysis/descriptive-analysis.md)

Compile to LaTeX with:

```
pandoc -s descriptive-analysis.md -o descriptive-analysis.tex
```

And add `\usepackage{morefloats}` to the LaTeX header and grep for `.png` to `.pdf`.

Compile to PDF with:

```
pdflatex descriptive-analysis
```
