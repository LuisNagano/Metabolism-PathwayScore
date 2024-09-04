# Metabolic Pathway Analysis in Lung Squamous Cell Carcinoma (LUSC)

## Overview
This project focuses on pathway scoring for Lung Squamous Cell Carcinoma (LUSC) using transcriptomic data. The analysis is based on the methodology from the paper ["Pan-cancer analysis of transcriptional metabolic dysregulation using The Cancer Genome Atlas"](https://doi.org/10.1038/s41467-018-07546-9) by Rosario et al., 2018, which was adapted and expanded for this analysis. The aim of this project is to identify key metabolic pathways that are dysregulated in LUSC samples compared to normal tissue, leveraging data downloaded from the [GDAC Firehose](https://gdac.broadinstitute.org).

## Objectives

1. Identify and quantify metabolic pathway dysregulation in LUSC compared to normal tissue.
2. Generate insights into the metabolic reprogramming occurring in LUSC.
3. Provide a foundation for potential therapeutic targets based on metabolic vulnerabilities.

## Data Source

- Raw count data: TCGA-LUSC dataset (https://gdac.broadinstitute.org)
- Initial dataset: TCGA-ACC raw counts file

## Methodology

1. **Data Preparation**:
   - Cleaned and processed raw count data from TCGA-LUSC.
   - Separated samples into Normal and Tumor groups.

2. **Differential Expression Analysis**:
   - Utilized limma package for differential gene expression analysis.
   - Calculated log fold changes and adjusted p-values.

3. **Pathway Scoring**:
   - Implemented a custom pathway scoring algorithm based on the PathwayScore method from Rosario et al. (2018).
   - Utilized 114 metabolic pathways from KEGG database.

   The pathway scoring process involves three crucial steps:

   a) **Gene Score Calculation**:
      The Gene Score is calculated using the following formula:

      ```
      Gene score = |log FC* - log(adj. p. val)|
      ```

      This score combines the magnitude of gene expression change (log fold change) with the statistical significance of that change (adjusted p-value). This approach allows us to prioritize genes that show both large and statistically significant changes in expression.

   b) **Pathway Score Calculation**:
      The Pathway Score is then calculated as:

      ```
      Pathway score = Σ(Gene scores) / √n
      ```

      Where n is the sample size for that particular tissue. This scoring method aggregates the impact of all genes within a pathway while accounting for differences in sample sizes across cancer types. This normalization is crucial for fair comparisons between different cancer cohorts.

   c) **Statistical Significance via Bootstrapping**:
      To determine the statistical significance of pathway scores, a bootstrapping approach is employed:
      - Pathway scores are randomly generated 100,000 times per pathway, based on the number of genes in the pathway.
      - The actual pathway score is then compared against this distribution of random scores.
      - A p-value is calculated based on where the actual score lies within the distribution.

      This bootstrapping method provides a robust way to assess the significance of pathway dysregulation, taking into account the size and structure of each pathway.

   The importance of this scoring system lies in its ability to:
   - Quantify the overall dysregulation of entire metabolic pathways.
   - Account for both the magnitude and statistical significance of gene expression changes.
   - Provide a standardized method for comparing pathway dysregulation across different cancer types and studies.
   - Identify significantly dysregulated pathways that may not be apparent from single-gene analyses.

4. **Statistical Analysis**:
   - Performed bootstrapping to determine statistical significance of pathway scores.
   - Adjusted pathway scores based on bootstrapping results.

5. **Visualization**:
   - Generated heatmaps and bar plots to visualize top dysregulated pathways.

## Key Findings

1. Identified the top 10 most dysregulated metabolic pathways in LUSC:
   - Nicotinamide Adenine Metabolism
   - Beta Alanine Metabolism
   - Phenylalanine Metabolism
   - Amino Sugar and Nucleotide Sugar Metabolism
   - Glycerolipid Metabolism
   - Oxidative Phosphorylation
   - Pyruvate Metabolism
   - Folate One Carbon Metabolism
   - Valine, Leucine and Isoleucine Degradation
   - Histidine Metabolism

![Top 10 LUSC Metabolic Pathways](https://github.com/LuisNagano/Metabolism-PathwayScore/blob/main/Figures/Top10_LUSC_Metabolic_Pathways.png)

![Top 10 LUSC Pathways (P-values)](https://github.com/LuisNagano/Metabolism-PathwayScore/blob/main/Figures/Top10_LUSC_Metabolic_Pathways_Pvalues.png)

2. Observed significant upregulation of the Pentose and Glucuronate Interconversions (PGI) pathway in LUSC, which is consistent with previous metabolomic studies showing increased UDP-D-glucuronate levels in lung cancer tissues.

3. Identified potential metabolic vulnerabilities that could be exploited for therapeutic interventions.

## Insights and Implications

1. The dysregulation of multiple amino acid metabolism pathways suggests altered protein synthesis and degradation processes in LUSC.

2. Upregulation of glycerolipid metabolism may indicate increased lipid synthesis to support rapid cell proliferation.

3. The significant dysregulation of oxidative phosphorylation could be related to the Warburg effect, a hallmark of cancer metabolism.

4. The alteration in folate one-carbon metabolism might affect nucleotide synthesis and DNA methylation processes, potentially contributing to genomic instability in LUSC.

## Future Directions

1. Validate key findings using metabolomics data from LUSC samples.
2. Investigate the potential of targeting highly dysregulated pathways for therapeutic interventions.
3. Expand the analysis to other cancer types to identify common and unique metabolic signatures.
4. Integrate this metabolic profiling with other omics data (e.g., proteomics, epigenomics) for a more comprehensive understanding of LUSC biology.

## Dependencies

- R (version 3.1.3 or higher)
- Bioconductor (version 3.1 or higher)
- R packages: limma, edgeR, ggplot2, Rsubread

## References

1. The Cancer Genome Atlas Research Network. (http://cancergenome.nih.gov/)
2. Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Research, 28(1), 27-30.
3. Rosario, S.R., et al. (2018). Pan-cancer analysis of transcriptional metabolic dysregulation using The Cancer Genome Atlas. Nature Communications, 9(1), 5330.

## Acknowledgments

This project is an adaptation of the PathwayScore algorithm described in the "Rosario et al., 2018" article. We acknowledge the original authors for their foundational work in metabolic pathway analysis.
