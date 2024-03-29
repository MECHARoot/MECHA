
V 2.1 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10893083.svg)](https://doi.org/10.5281/zenodo.10893083)

V [2.0](https://github.com/MECHARoot/MECHA/releases/tag/2.0)

V [1.0](https://github.com/MECHARoot/MECHA/releases/tag/1.0)

# MECHA

## Description

MECHA is an explicit cross-section model of the root hydraulic anatomy which connects hydraulic concepts across scales.

The model computes horizontal water flow at the level of individual cells, quantifies the contribution of water composite pathways, and predicts root radial permeability (kr), using detailed anatomical descriptions and experimental data on the permeability of cell walls (kw), membranes (Lp) and the conductance of individual plasmodesmata (KPD).

## Installation

Conda/Mamba

```{bash}
conda env create -f environment.yml

conda activate mecha_env
```

## Use

MECHA uses four parameters files as input, stored in the /in folder in the Projects MECHA directory:

Projects/<name_of_the_project>/in/

General.xml : General simulation parameters
BC.xml : Boundary conditions for the simulation. This is where you choose wether the soil is dry, wet, or only partially touch the root.
Geometry.xml : Parameters for the root section geometry. This is where you define which root section to use in the simulation.
Hydraulics.xml : Hydraulic parameters of the different variables in the simulation

#### Run MECHA
```{bash}
python MECHA.py
```

OUPUT FILES
MECHA output files in a folder defined in the General.xml input file, under the Output tag. Different type of files are generated:

*.PVTK files : These are geometry files that can be open in Paraview
Macro_prop_i.txt files : Theses files contain synthetic informations about the simulation results. They contain the radial data.

## Cite

Valentin Couvreur, Marc Faget, Guillaume Lobet, Mathieu Javaux, François Chaumont, Xavier Draye, Going with the Flow: Multiscale Insights into the Composite Nature of Water Transport in Roots, Plant Physiology, Volume 178, Issue 4, December 2018, Pages 1689–1703, https://doi.org/10.1104/pp.18.01006

## Licence

MECHA is released under a GPL-3 licence, which means that redistribution and use in source and binary forms, with or without modification, are permitted under the GNU General Public License v3 and provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
