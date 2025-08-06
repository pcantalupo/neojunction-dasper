# Neojunction Nextflow Pipeline
This Nextflow pipeline performs neojunction analysis on individual samples by using [dasper](https://github.com/dzhang32/dasper) to annotate RNA-seq splice junctions. It processes STAR splice junction files, applies a filtering step, and combines the filtered junctions for neojunction analysis. The pipeline generates a report with the proportion of known versus novel junctions per sample, along with a percentage breakdown of each type of neojunction.

---
## Requirements

- [Nextflow](https://www.nextflow.io/) (developed with 24.10.5)
- R with required packages installed
    - Important: the dasper package has has been removed from Bioconductor and no longer installable with current R or Bioconductor releases.
    - To install dasper, you need an old version of R (tested with 4.0.0) and Bioconductor (tested with 3.12)
- Container engine (Recommended) such as Docker or Singularity installed
- STAR-generated SJ.out.tab files for each sample
- GTF annotation file for your reference genome

---

## Input Files

- **Samplesheet (`samplesheet.csv`)**: CSV file with headers containing:
  - `sample_id`: Unique sample name
  - `sj_file`: Path to the STAR splice junctions file (SJ.out.tab)
  - `library_size`: Library size value for scaling junction counts

Example:

```csv
sample_id,sj_file,library_size
S1,star/sample1_SJ.out.tab,1000000
S2,star/sample2_SJ.out.tab,2000000
```

- **GTF file**: Reference annotation file corresponding to your genome build.

---

## Usage
Basic command to run the pipeline:
```bash
nextflow run pcantalupo/neojunction-dasper --samplesheet samplesheet.csv --gtf path/to/annotations.gtf
```

To run the pipeline with a container engine. Singularity and docker are supported.
```bash
nextflow run pcantalupo/neojunction-dasper -profile <singularity|docker> --samplesheet samplesheet.csv --gtf path/to/annotations.gtf
```

---

## Output

- Results are saved in `results/` by default (`--outdir`).
- Dasper output and dasper filtered files will be in `results/dasperannots` by default (`--dasper_outdir`).
- Neo-junction results will be saved as `results/neojunction_results.tsv`.
- Nextflow pipeline information is available in the pipeline info directory (`results/pipeline_info`).

---
## Testing
See the `tests/` directory for a test dataset to check your installation.

