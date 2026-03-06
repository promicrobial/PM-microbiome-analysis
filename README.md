# Microbiome Analysis of Iron/Zinc-Biofortified Pearl Millet Trial

This repository contains the analysis pipeline for a randomized controlled feeding trial investigating the effects of iron/zinc-biofortified pearl millet (FeZnPM) on the gut microbiome of children in Mumbai, India (Clinical Trials.gov ID: [NCT02233764](https://clinicaltrials.gov/study/NCT02233764?term=NCT02233764&rank=1)).

## Repository Structure

```
├── analysis/     # Analysis Quarto notebooks organized by topic
├── _assets/      # Style files and images
├── manuscript/   # Manuscript files and references
├── scripts/      # Analysis and processing scripts
└── utils/        # Utility functions and helpers
```

## Analysis overview

Analysis was conducted in R 4.4.2 and package requirements and versions are available in the renv.lock file.

- Quality control and preprocessing of metagenomic data
- Classification by alignment to WoL2 reference phylogeny
- Diversity analyses (alpha and beta)
- Differential abundance testing
- Functional analysis (KEGG pathways and GO terms)
- Iron status subgroup analysis

## Citation

Huey, S.L., Cole, N.L., Pagani, I., González, A., Finkelstein, J.L., Haas, J.D., Udipi, S.A., Ghugre, P., Potdar, R.D., Knight, R., Mehta, S., (in preperation). Effect of a complementary feeding intervention on the gut microbiota in 12–18-month-old children.
