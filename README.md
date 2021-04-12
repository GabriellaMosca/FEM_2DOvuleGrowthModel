# FEM_2DOvuleGrowthModel

[![DOI](https://zenodo.org/badge/233677761.svg)](https://zenodo.org/badge/latestdoi/233677761)

2D continuous models of ovule growth based on FEM

To run the models you need to have MorphoMechanX installed. To do so, go to www.morphomechanx.org (select again the MorphoMechanX tab) and follow the instructions provided there. As a reference, this model is guaranteed to run against version 2.0.2-858 (git commit: 8eb2aec9055221e3e351ffe064bd059186717374)

Once MorphoMechanX is installed, download this repository to your local computer. Open a terminal and move into this repository folder ("yourPath"+ FEM_2DOvuleGrowthModel).
From there type: 
- make clean
- make run

This will launch MorphoMechanX with the model already loaded (you can see the model Processes under the Tab: Process/Model). 

To load and run a specific model:
- enter in one of the subfolders of this repository (e.g. FEM-Model1);
- find the MorphoDynamX view file, i.e. the file terminating with ".mdxv" (called "FemMembranes.mdxv") and drag and drop it into the MorphoMechanX window
- within the MorphoMechanX GUI, go to the Process/Model tab
- double click on "Model/CCF/01 FEM Membranes"

In this way the specific growth simulation will start. Mesh templates will be saved at the simulation times indicated in the field "Simulation time to save mesh" under the Process: "Model/CCF/01 FEM Membranes".

For further details on how to edit model parameters or their meaning, refer to the Appendix 1 provided with the paper (Hernandez-Lagana, Mosca, Mendocilla-Sato et al.). 
The Appendix can be found as pre-print also here https://www.biorxiv.org/content/10.1101/2020.07.30.226670v3.supplementary-material under "Supplemental File 1"

For further questions, contact: gabriella.mosca@gmail.com
