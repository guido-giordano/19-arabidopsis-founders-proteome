# Arabidopsis Founders Proteome Pipeline

Reproducible Python pipeline for the analysis of the 19 Arabidopsis thaliana founder lines proteome.

## Overview
This project implements a structured workflow to process, integrate, and analyse quantitative proteomics data generated from multiple Arabidopsis accessions. Developed in at the University of Cambridge (Lilley Lab), in the context of the project "What determined protein stability in plants", focusing on understanding the determinants of protein stability in Arabidopsis thaliana

## Features
- Quality control of large-scale proteomics datasets (filtering, consistency checks, summary metrics)
- Construction of orthology-informed protein matrices across multiple accessions
- Integration of quantitative protein abundance data (iBAQ-based)
- Automated generation of publication-ready visualisations
- Handling of multi-file datasets and large matrices typical of modern proteomics workflows

## Technical Highlights
- Developed in Python (pandas, numpy, matplotlib, seaborn)
- Designed for reproducible and modular data analysis
- Compatible with remote execution on Linux/HPC environments
- Structured to support scalability to multi-omics datasets

## Future Extensions
- Integration with transcriptomics datasets
