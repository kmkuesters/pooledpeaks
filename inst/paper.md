---
title: 'pooledpeaks: Open-Source Fragment Scoring and Analysis Workflow for Population 
Genetic Analyses of Microsatellite Markers in Pooled Samples'

tags:
  - R
  - population genetics
  - fragment analysis
  - pooled samples
  - .fsa
authors:
  - name: Kathleen M. Kuesters
    orcid: 0000-0003-0238-7889
    corresponding: true 
    affiliation: 1 
  - name: Jessica M. Blanton
    orcid: 0000-0002-1653-7836
    affiliation: 1
  - name: Jeffrey D. Kovach
    orcid: 0000-0001-9051-7676
    affiliation: "2,3"
  - name: Walter A. Blank
    affiliation: 4
  - name: Lúcio M. Barbosa
    orcid: 0000-0001-7272-0515  
    affiliation: "5,6"
  - name: Luciano K. Silva
    orcid: 0000-0002-7739-4082
    affiliation: 6
  - name: Mitermayer G. Reis
    orcid: 0000-0002-3051-9060
    affiliation: "6,7,8"
  - name: Ronald E. Blanton
    orcid: 0000-0001-6655-7336
    affiliation: 1
affiliations:
 - name: Department of Tropical Medicine and Infectious Disease, Tulane University Celia Scott Weatherhead School of Public Health and Tropical Medicine, New Orleans, LA 70112, United States of America
   index: 1
   ror: 04vmvtb21
 - name: Case Western Reserve University School of Medicine, Cleveland, OH 44106, United States of America
   index: 2
   ror: 051fd9666
 - name: Department of Quantitative Health Sciences, Lerner Research Institute, Cleveland Clinic, Cleveland, OH, United States of America
   index: 3
   ror: 03xjacd83
 - name: Independent Researcher, United States of America
   index: 4
 - name: Bahiana School of Medicine and Public Health, Av. Silveira Martins, n 3386, Salvador, Bahia, 41150-100, Brazil
   index: 5
   ror: 0300yd604
 - name: Instituto Gonçalo Moniz, Fundação Oswaldo Cruz, Ministério da Saúde, Rua Waldemar Falcão, 121, Candeal, CEP 40296-710, Salvador, Bahia, Brazil
   index: 6
 - name: Faculdade de Medicina, Universidade Federal da Bahia, Praça XV de novembro, s/n - Largo do Terreiro de Jesus, CEP 40026-010, Salvador, Bahia, Brazil
   index: 7
   ror: 03k3p7647
 - name: Department of Epidemiology of Microbial Diseases, Yale School of Public Health, 60 College St, New Haven, Connecticut, 06510, United States of America
   index: 8
   ror: 03v76x132
date: 1 February 2025
bibliography: pooledpeaks_paper.bib
---

# Summary

Microsatellite markers are short, highly variable, multi-repeat DNA sequences 
(aka short tandem repeats) that appear throughout the genome and can be used to
estimate population genetic metrics [@silva; @vieira]. 
These markers are frequently evaluated using fragment analysis which is based 
on Sanger sequencing. The `pooledpeaks` R package provides tools to analyze 
fragment analysis results (.fsa files). It provides functions that fall in 
three subcategories: 1) peak scoring, 2) data manipulation, and 3) genetic 
analysis. The package was designed for the use of microsatellite markers on 
pooled parasite samples, but the peak scoring functions are applicable to any 
fragment analysis. The peak scoring functions were partially adapted from 
Fragman, a package designed to score microsatellite markers in cranberries
[@covarrubias]. Although Fragman works for the older file version, newer
versions cannot be read. In addition to revising this outdated function, we 
also added features including expanded scoring parameter options and exporting 
resulting scoring plots as a pdf file for review. The data manipulation 
functions were created to clean and format the data from the called peaks and 
transform them into allele frequencies. These frequencies can then be input 
into the genetic analysis functions for calculation of diversity and 
differentiation measures adapted from a range of papers
[@long;@jost;@nei;@foulley;@chao]. An in-depth 
walk-through of how to use the analysis pipeline can be found in the vignette.

# Statement of Need

While a plethora of methods exist for downstream statistical analysis of allele
frequencies, processing raw fragment data is limited by available software. Of 
the limited software that can read the .fsa binary raw data file format, nearly
all require purchase or registration, are primarily built for Windows, are 
inefficient for analyzing large batches of files, and are highly dependent on 
individual researcher experience. Additionally, a previous R package allowing 
for the analysis of .fsa files is incompatible with the updated file version. 
When using fragment analysis for microsatellite markers on pooled samples, once
the raw data is extracted and scored, it must be cleaned and transformed into 
allele frequencies using a second software, such as Excel, which is limited in
its capacity for automation and version control. Another platform shift is 
often required to analyze the resulting allele frequencies. These factors 
highlight the need for a comprehensive scoring and analysis pipeline that is 
open-source, offline, reproducible, consistent between researchers, and that 
does not require platform switching between steps. 

# Ongoing Research Projects 
This package is currently being used to analyze genetic clustering of 
*Schistosoma mansoni* pooled egg samples from four Brazilian communities, as 
well as the relatedness of *Schistosoma haematobium* populations around Lac de 
Guiers in Senegal and from Gabon.

# Acknowledgements (Financial)

This work was financially supported by the NIH as part of 1R01AI121330.

# References
