# Evidence-of-widespread-sympatry-as-a-diversification-driver-for-Emberizoidea-Aves-Passeriformes-
This repository contains the necessary tools and data to allow the recreation of the analyses performed in this study.


File description:
Data: This folder contains the core files used for this study
	- BAMM-MMC_taxonomy.csv: Clements et al., (2019) taxonomic classification for all the species used in this study
	- breeding_phylip_k8.txt: presence and absence file in PHYLIP format, obtained using the phyloregion method, used for the ancestral area reconstruction of the Emberizoidea superfamily.
	- Breeding-phylogeny_MCC.tre: Barker et al., (2015) Maximum clade credibility phylogeny.
	- Embemaps_v1.Rdata: Spatial Polygon DataFrame of the Emberizoidea species' global distribution
	- family phy.tre: Phylogenetic tree with a single tip per family, used for the PGLS models
	- nodedtree.txt: Barker et al., (2015) Maximum clade credibility phylogeny with named nodes.
	- parprobx.csv: nodedtree transformed into a Data.Frame, used for the transition rate calculations
	- pglsread_rage: Data Frame with the data per family used for the PGLSs, including mean speciation rate (BAMM), Local Extinction rate (eTR; extinction_rate), Dispersal rate (dTR; dispersal_rate), family age (age), total richness per family (total_richness).

Scripts: This folder contains the annotated scripts necessary to perform the analyses used in this study.
	- masterscript_phyloregion.R: This script contains the code required to calculate the phylogenetic regions for Emberizoidea using the phyloregion package for R.
	- masterscript_biogeobears-rates.R: This script contains the code required to reconstruct the ancestral area models, calculate the Stochastic Biogeographic Mappings (BSM), identify diversification shifts and calculate the Transition and Speciation rates (Including lambdaBAMM, DR, ClaDS and Method of Moments) used in this study.
	- masterscript_pgls.R: This script contains the instructions to perform the PGLSs in this study.

phyloregion: This folder contains the phylogenetic distance matrix for Emberizoidea to save time; calculated in mastrescript_phyloregion and loaded when runslow is set to F. Otherwise serves as a container for all the phyloregion results.
	- phylobeta.Rdata: Simpson's phylogenetic index distance matrix.

biogeobears: This folder contains the biogeographic models performed by BioGeoBEARS, using the phylogenetic regions calculated using K8, also serves as a container for all the BioGeoBEARS results.
	- Allrange_DEC_k8_unconstrained_v2.Rdata: DEC model using 8 distinct regions for Emberizoidea
	- all range_DEC+J_k8_unconstrained_v2.Rdata: DEC model considering jump dispersal
	- All range_DIVALIKE_k8_unconstrained_v2.Rdata: DIVA model using 8 distance regions for Emberizoidea
	- All range_DIVALIKE+J_k8_unconstrained_v2.Rdata: DIVA model considering jump dispersal
	- All range_BAYAREALIKE_k8_unconstrained_v2.Rdata: BayArea model using 8 distance regions for Emberizoidea, used for the rest of the analyses.
	- All range_BAYAREALIKE+J_k8_unconstrained_v1.Rdata: BayArea model considering jump dispersal

bsm: This folder contains a single object with the basic inputs for the calculation of the stochastic biogeographic mapping (SBM), using the BayArea model and can be used to contain the SBM results. The inputs object can be created with the masterscript_biogeobears-rates.R script.
	- BSM_inputs_file.Rdata: Inputs necessary to run the SBMs, calculated for the BayArea model

Diversification: This folder contains the posterior distributions of the speciation rate calculations, and the missing species for each family in the Barkers et al. (2015) phylogeny according to the Clemens Taxonomy.
	- BAMM: data used for the BAMM calculations
		- BAMM_chain_swap.txt: frequency of hot/cold chain swapping and chain convergence
		- BAMM_event_emberizo_1.txt: Posterior distributions of all parameters associated with the macroevolutionary rate regimes
		- BAMM_mcmc_emberizo_1.txt: file containing basic attributes of the MCMC chain 
		- BAMM_run_emberizo_1.txt: Summary of the BAMM analysis.
	- ClaDS: Posterior distributions of ClaDS
		- clads_output.RData: All the parameters associated with the macroevolutionary rate regimes rate regimes, calculated using the ClaDS method.
	- rates: This folder contains the information of the missing species per family in the Emberizoidea phylogeny. This information is required for the calculation of the Method of Moments diversification rate
		- missing_families.csv: information of the missing species per family.
