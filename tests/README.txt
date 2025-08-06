Testing the Pipeline Locally
This directory contains a test dataset (2 samples, chr22 only) to verify that the pipeline runs correctly

Running the Test
1. Run the pipeline locally
Make sure all required software dependencies are installed on your machine. Then execute:
nextflow run ../main.nf --samplesheet samplesheet.csv --gtf data/chr22.gtf

2. Run the pipeline using a container
For Singularity:
nextflow run ../main.nf -profile singularity --samplesheet samplesheet.csv --gtf data/chr22.gtf

For Docker:
nextflow run ../main.nf -profile docker --samplesheet samplesheet.csv --gtf data/chr22.gtf

Verifying Results
- After successful execution, output files are in the results/ directory
- Check that the generated output matches the expected results by comparing md5 checksums

md5sum expected_results.tsv results/neojunction_results.tsv
eccc5881838c6b1b3051878787470409  expected_results.tsv
eccc5881838c6b1b3051878787470409  results/neojunction_results.tsv

Troubleshooting
- If you get file not found errors, double-check your working directory and the paths in the samplesheet.
- When using containers, ensure you have network access and the correct profile is activated.
- Consult .nextflow.log and the work directory under work/ for detailed process logs if the pipeline fails.

