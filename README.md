*S. faciatus* RNA-seq analysis

**Main author:** Christelle Leung  
**Affiliation:** Fisheries and Oceans Canada (DFO)  
**Group:** Laboratory of genomics  
**Location:** Maurice Lamontagne Institute  
**Affiliated publication:** Leung C, Guitard J, Senay C, Bourret A,
Chabot D and Parent GJ (In prep). Past environments modulate response to
fluctuating temperatures in a marine fish species.  
**Contact:** e-mail: <christelle.leung@uqtr.ca>

- [Objective](#objective)

- [Summary](#summary)

- [Status](#status)

- [Contents](#contents)

- [Requirements](#requirements)

- [Acknowledgements](#acknowledgements)

## Objective

Assessing *Sebastes fasciatus* temperature-related transcriptional
plasticity

## Summary

Predicted ocean warming will impact the survival and structure of
various marine organisms, in particular ectotherms. Phenotypic
plasticity enables species to cope with environmental changes, providing
a vital buffer for evolutionary adaptations. Yet, the dynamics and
molecular mechanisms underpinning these plastic responses remain largely
unexplored. Here, we assessed the impact of temperature acclimation on
organisms’ capacity for thermal plasticity. We conducted a genome-wide
transcriptomic analysis on the Acadian redfish, *Sebastes fasciatus*,
exposed to four temperatures (2.5, 5.0, 7.5 and 10.0 ℃) over a long-term
period (up to 10 months) followed on some individuals by a short-term
temperature change (+2.5 °C or -2.5 °C for 24 hours), simulating natural
temperature variation the species could encounter. Our results showed a
dynamic transcriptional response to temperature involving various gene
functions. The rapid response to temperature shifts, coupled with the
sustained expression of specific genes over an extended period
highlighted the species’ capacity for plasticity in response to
temperature changes. We also detected a significant effect of the
interaction between the long- and short-term temperature exposures on
gene expression, highlighting the influence of the past environment on
the response to short-term temperature changes. Specifically, fish
acclimated to higher temperatures demonstrated an increased
stress-related response to environmental fluctuations, as evidenced by
both the shape of their reaction norms and the involvement of
stress-related gene functions. This result suggests that temperature
conditions predicted for the near future in the Northwest Atlantic will
trigger reduced adaptive plasticity to environmental fluctuations,
highlighting the species’ vulnerability to ocean warming.

## Status

Ongoing

## Contents

00_Data: Folder to store raw data, metaData and intermediates files -
Large files are stored locally.

01_Scritps: All scripts used to perform the entire analysis,
including:  
- 01_RAW_quality_assessment.R: Check for raw read quality (FastQC &
MultiQC)  
- 02_Trim_reads.R: Trim reads for low quality and adaptors (Trim
Galore!)  
- 03_Mapping.R: Reads mapping to the reference genome (Hisat2)  
- 04_Assembling: *De novo* transcritpome construction, guided by
reference genome (StringTie) - Structural genome annotation  
- 05_Counts.R: Get the gene/transcript counts (FeatureCounts)  
- 06_FunctionAnnot.R: Functional annotation of *de novo* transcritpome  
- 07_DETs_stats.R: Variation partioning and differential gene expression
analyses  
- 08_Sex_effect.R: Test for sex effect on gene expression  
- 99_Search_GOterms.R: Script to search specific GO terms in
differential gene expression results

02_Results: Intermediate and final results: *de novo* transcritpome,
count and statistical result tables, figures - Large files are stored
locally.

## Requirements

Raw sequence data for the *S. fasciatus* genome assembly is available
under the Genbank accession \#JBJQUQ000000000, the RNA-seq and ddRAD-seq
datasets are available in the Sequence Read Archive (SRA) under the
BioProject accession \#PRJNA1208449

## Acknowledgements

We thank D.C. Chavarria, J. Gagnon, D. Picard, T. Hansen, F. Hartog, K.
MacGregor, JD. Tourangeau-Larivière and J. Heinerth for technical
support during the experiment at the DFO Maurice Lamontagne Institute
(Mont-Joli, QC, Canada) and for fish collection. We also acknowledge G.
Bardaxoglou, J. Larivière and G. Cortial for DNA and RNA extraction.
