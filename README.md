# An evaluation workflow of Orthanq for virus applications.

## Configuration 

The `config.yaml` of the main workflow should be configured with:
-  single sample is simulated, `simulate_given: true and simulate_pandemics: false` and `lineages.tsv` should contain lineages and corresponding fractions with lineage names provided with file names under `resources/lineages`.
-  pandemics simulation, `simulate_given: false and simulate_pandemics: true` and `number_of_samples` should be set.
-  analysis of given sra IDs: `simulate_given: false and simulate_pandemics: false` and `samples_wastewater.tsv` should be filled in.

In addition, under `config/pep`, `config_filled.yaml` should be filled with the correct sample sheet depending on a simulation of a single sample, pandemics simulation or analysis of given SRA IDs. Sample sheets should be filled in parallel with the mode of workflow execution, under `config/pep` as `samples_filled_single_sample.tsv`, `samples_filled_pandemics.csv` and `samples_filled_sra.csv`.