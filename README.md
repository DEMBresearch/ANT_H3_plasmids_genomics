# ANT_H3_plasmids_genomics
This repo includes data and code used for the comparative genomic analysis of *Psychrobacter* sp. ANT_H3 plasmids compared with other plasmids from *Psychrobacter* spp. described in the article "Structure and functions of a multireplicon genome of Antarctic *Psychrobacter* sp. ANT_H3 – characterization of the genetic modules suitable for the construction of the plasmid-vectors for cold-active bacteria".

## How to run the analysis
1. Create conda environment with all required software

You can either use the provided YAML file or install all the required tools by yourself
```bash
# using conda
conda env create --file anth3_env.yml
# or mamba which resolves dependencies more efficiently
mamba env create --file anth3_env.yml
```
2. Activate the conda env
```bash
conda activate anth3
```

3. Run the BASH script where all the commands have comments
```bash
bash run_analysis.sh
```

4. All output files will be in `analysis` directory.
- `analysis/circos/ANT_H3.nucl.svg` file was further processed to create panel A in Figure 2.
- `analysis/network/net_noG_simple.graphml` file was opened with [Gephi](https://github.com/gephi/gephi/) v0.9.7 and further filtered and processed manually to generate panel B in Figure 2.