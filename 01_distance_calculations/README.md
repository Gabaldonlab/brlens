# Distance calculation module

This module allows the user to calculate the tip-to-internode and
tip-to-tip distances.

## Structure of the module

```
01_distance_calculations/
├── get_t2in_dist.py: tip-to-internode distances calculation
├── get_t2t_dist.py: tip-to-tip distances calculation
├── treefuns.py: functions to work with trees
└── utils.py: general functions to manage files and folders
```

## Usage
### Data preparation
To run a distances calculation analysis the first step is to develop
a dataset of trees and a table that relates species to events. The main
files the user will need will be:

1. `node_data.tsv`: this file contains the species code (common in the trees and tables) belonging to the groups of interests (can be numbers as in the following example or names) in columns:

```
26	  27	30	44
NA	  NA	  NA	  NA
MONDO	NA	  NA	  NA
CHOHO	CHOHO	NA	  NA
ECHTE	ECHTE	NA	  NA
LOXAF	LOXAF	NA	  NA
ERIEU	ERIEU	ERIEU	NA
PTEVA	PTEVA	PTEVA	NA
FELCA	FELCA	FELCA	NA
PANPR	PANPR	PANPR	NA
CANFA	CANFA	CANFA	NA
HORSE	HORSE	HORSE	NA
PIG	  PIG	  PIG	  NA
SHEEP	SHEEP	SHEEP	NA
TUPGB	TUPGB	TUPGB	NA
RABIT	RABIT	RABIT	NA
SPETR	SPETR	SPETR	NA
JACJA	JACJA	JACJA	NA
MOUSE	MOUSE	MOUSE	NA
RAT	  RAT	  RAT	  NA
PROCO	PROCO	PROCO	PROCO
CARSF	CARSF	CARSF	CARSF
PAPAN	PAPAN	PAPAN	PAPAN
MACMU	MACMU	MACMU	MACMU
HUMAN	HUMAN	HUMAN	HUMAN
```

2. `sp_tree.nwk`: the species tree using the species codes in newick format.

```
(ORNAN:1.997583,(MONDO:1.421804,((CHOHO:0.78658,(ECHTE:0.658754,LOXAF:0.658755):0.127826):0.014895,((ERIEU:0.697734,(PTEVA:0.6811,(((FELCA:0.154884,PANPR:0.154885):0.364435,CANFA:0.51932):0.15508,(HORSE:0.662792,(PIG:0.549832,SHEEP:0.549831):0.112961):0.011608):0.0067):0.016634):0.053029,((TUPGB:0.689588,(RABIT:0.661998,(SPETR:0.60782,(JACJA:0.459977,(MOUSE:0.104451,RAT:0.104451):0.355526):0.147844):0.054178):0.02759):0.008291,(PROCO:0.648951,(CARSF:0.615783,((PAPAN:0.104655,MACMU:0.104654):0.199438,HUMAN:0.304093):0.31169):0.033168):0.048929):0.052883):0.050712):0.620329):0.575779);
```

3. `trees.nwk`: the gene trees are stored in a tabular format containing the seed sequence, the model used for its infence, the log-likelihood of the tree and the newick formatted tree:

```
Phy0007XW0_HUMAN	JTTDCMut+R4	-14213.7528	(Phy003IMAD_MACMU:0.0000000000,(Phy00FD1NC_PAPAN:0.0000024491,...
```

### Running distances calculation
#### Tip-to-internode distances
The script arguments can be obtained by this `python3 get_t2in_dist.py -h`:
```
Usage: get_t2in_dist.py [options]

Options:
  -h, --help            show this help message and exit
  -i <file.nwk>, --input=<file.nwk>
                        File with multiple trees in newick, each line has the
                        format: seed      model   likelihood      newick
  -s <file.nwk>, --sptree=<file.nwk>
                        Newick file containing the species tree.
  -c <file.tsv>, --clades=<file.tsv>
                        Species belonging to clades, columns show the clades
                        and rows the species belonging to them.
  -l <SPECIES>, --seedsp=<SPECIES>
                        Species' code for the seed.
  -n <group_column_name>, --normgroup=<group_column_name>
                        Normalisation group header in clades dataframe.
  -t <N>, --threads=<N>
                        Number of threads.
  -p </path/to/dir/prefix> or <prefix>, --prefix=</path/to/dir/prefix> or <prefix>
                        Output prefix.
  -r, --redo            Ommit done file and redo.
```

To run the command using the files described before:
```bash
python3 get_t2in_dist.py -i trees.nwk -s sp_tree.nwk -c node_data.tsv -l HUMAN -n 44 -p outputs -t 2
```

#### Tip-to-tip distances
The script arguments can be obtained by this `python3 get_t2in_dist.py -h`:
```
Usage: get_t2t_dist.py [options]

Options:
  -h, --help            show this help message and exit
  -i <file.nwk>, --input=<file.nwk>
                        File with multiple trees in newick, each line has the
                        format: seed      model   likelihood      newick
  -c <file.tsv>, --clades=<file.tsv>
                        Species belonging to clades, columns show the clades
                        and rows the species belonging to them.
  -n <group_column_name>, --normgroup=<group_column_name>
                        Normalisation group header in clades dataframe.
  -t <N>, --threads=<N>
                        Number of threads.
  -p </path/to/dir/prefix> or <prefix>, --prefix=</path/to/dir/prefix> or <prefix>
                        Output prefix.
  -r, --redo            Ommit done file and redo.
```

To run the script on your `trees.nwk` file:
```bash
python3 get_t2t_dist.py -i trees.nwk -c node_data.tsv -n 44 -t 2
```

For obtaining a filtered set of `R` `data.frames` which can be used in
the inference process, you have to run the [distances_filtering.R](distances_filtering.R):

```bash
Rscript distances_filtering.R
```

It creates a set of `RData` files which can be imported in the inference
module.