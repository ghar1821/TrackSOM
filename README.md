# TrackSOM

TrackSOM is a software protocol for creating maps of how the dynamic immune response evolves over time and/or disease severity.
These data are captured in a sequence of cytometry datasets representing points in time, and/or disease stages.
TrackSOM will identify the immune cell populations (aka phenotypes) in these data, and how they evolve and vary over the dataset sequence.
This includes branching/splitting from common progenitors, appearing and disappearing, changing in absolute or relative phenotypic cell counts, or movements through marker space that signify differentiation and/or functional changes.
TrackSOM includes a novel clustering and cluster-tracking algorithm, and a variety of visualisation techniques to enable interpretation of the resultant maps.

## Step by step walkthrough

A step by step R markdown notebook is available in the [`inst/sample_workflow/TrackSOM_workflow.md`](https://github.com/ghar1821/TrackSOM/blob/master/inst/sample_workflow/TrackSOM_workflow.md) file.

A complete R script with the workflow is also available in the [`inst/sample_workflow/TrackSOM_workflow.R`](https://github.com/ghar1821/TrackSOM/blob/master/inst/sample_workflow/TrackSOM_workflow.R) file.

Both sample workflows make use of the synthetic dataset reported in our manuscript.

## System Requirements
- Hardware requirements: Standard computer. TrackSOM has been tested on a Macbook with >= 8.0GB RAM.
- OS requirements: MacOS or Windows. TrackSOM has been tested on MacOS Big Sur 11.4 and Windows 10.
- R package dependencies (TrackSOM has been tested on the dependencies' version listed below):
	- R (3.6.0)
	- pdist (1.2)
	- rlist (0.4.6.1)
	- FlowSOM (1.18.0)
	- flowCore (1.52.1)
	- data.table (1.13.6)
	- gtools (3.9.2)
	- tidygraph (1.2.0)
	- ggraph (2.0.5)
	- RColorBrewer (1.1-2)
	- viridis (0.6.1)
	- stringr (1.4.0)
	- igraph (1.2.6)
	- ggplot2 (3.3.5)
	- Biobase (2.46.0)

## Installation guide
Installation is currently only available through the `devtools` R package.
```
install.packages("devtools")
devtools::install_github("ghar1821/TrackSOM")
```

TrackSOM will be made available in Bioconductor soon.

Typical install time < 1 minute.

This `TrackSOM` repository captures the software protocol only, to reduce the download size for users.
The datasets used in our manuscript are [stored at the Open Science Framework](https://osf.io/8dvzu/) (they are too large for GitHub), and the [corresponding evaluation scripts use in the manuscript are available here](https://github.com/ghar1821/TrackSOM-evaluations).

## Manuscript

R scripts to reproduce evaluations and figures in our paper introducing the `TrackSOM` algorithm are available on [`TrackSOM-evaluations`](https://github.com/ghar1821/TrackSOM-evaluations) repository.

Please cite the following preprint if you find TrackSOM useful in your research.

> Givanna H. Putri, Jonathan Chung, Davis N. Edwards, Felix Marsh-Wakefield, Suat Dervish, Irena Koprinska, Nicholas J.C. King, Thomas M. Ashhurst and Mark N. Read. (2021). TrackSOM: mapping immune response dynamics through sequential clustering of time- and disease-course single-cell cytometry data. bioRxiv 2021.06.08.447468; doi: https://doi.org/10.1101/2021.06.08.447468

The preprint titled **TrackSOM: mapping immune response dynamics through sequential clustering of time- and disease-course single-cell cytometry data** can be found on [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.08.447468v1).

The citation in bibtex format:
```
@article {Putri2021.06.08.447468,
	author = {Putri, Givanna Haryono and Chung, Jonathan and Edwards, Davis N and Marsh-Wakefield, Felix and Dervish, Suat and Koprinska, Irena and King, Nicholas JC and Ashhurst, Thomas Myles and Read, Mark Norman},
	title = {TrackSOM: mapping immune response dynamics through sequential clustering of time- and disease-course single-cell cytometry data},
	elocation-id = {2021.06.08.447468},
	year = {2021},
	doi = {10.1101/2021.06.08.447468},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Mapping the dynamics of immune cell populations over time or disease-course is key to understanding immunopathogenesis and devising putative interventions. We present TrackSOM, an algorithm which delineates cellular populations and tracks their development over a time- or disease-course of cytometry datasets. We demonstrate TrackSOM-enabled elucidation of the immune response to West Nile Virus infection in mice, uncovering heterogeneous sub-populations of immune cells and relating their functional evolution to disease severity. TrackSOM is easy to use, encompasses few parameters, is quick to execute, and enables an integrative and dynamic overview of the immune system kinetics that underlie disease progression and/or resolution.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2021/06/09/2021.06.08.447468},
	eprint = {https://www.biorxiv.org/content/early/2021/06/09/2021.06.08.447468.full.pdf},
	journal = {bioRxiv}
}
```

## Support and contribute
We are continuously working to improve TrackSOM and welcome feedbacks.
If you have questions or have encountered issues using TrackSOM, please post an [`issue ticket`](https://github.com/ghar1821/TrackSOM/issues).

## Publications using TrackSOM
TrackSOM was used to analyse COVID-19 temporal data published in Cell Reports Medicine journal:


> Koutsakos, M., Rowntree, L.C., Hensen, L., Chua, B.Y., van de Sandt, C.E., Habel, J.R., Zhang, W., Jia, X., Kedzierski, L., Ashhurst, T.M. and Putri, G.H., 2021. Integrated immune dynamics define correlates of COVID-19 severity and antibody responses. Cell Reports Medicine, p.100208.

R scripts to reproduce evaluations and figures in the paper are available on [`TrackSOM-evaluations`](https://github.com/ghar1821/TrackSOM-evaluations) repository.

## License
TrackSOM is released under the GPL-3.0 lincense, included [here](LICENSE).
