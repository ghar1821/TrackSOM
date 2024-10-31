# TrackSOM

TrackSOM is an algorithm for investigating how the dynamic immune response evolves over time and/or disease severity.
It maps out the temporal development of immune response captured in a sequence of cytometry datasets representing points in time or disease stages.
TrackSOM will identify the immune cell populations (aka phenotypes) in these data, and track how they evolve over the dataset sequence.
This includes branching/splitting from common progenitors, appearing and disappearing, changing in absolute or relative phenotypic cell counts, or movements through marker space that signify differentiation and/or functional changes.
TrackSOM includes a novel clustering and cluster-tracking algorithm, and a variety of visualisation techniques to enable interpretation of the resultant maps.

## Step by step walkthrough

A sample workflow is available on the [ImmuneDynamics website](https://immunedynamics.io/spectre/tutorials/time_series/TrackSOM-workflow.html) or from the Articles tab on [TrackSOM website](https://ghar1821.github.io/TrackSOM/).

The raw sample workflow files are also available in the [`doc`](https://github.com/ghar1821/TrackSOM/blob/master/doc).
They are available in `.Rmd`, `.R`, and `.html` format.

Approximate run time for the workflow: 1-2 minutes.

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

Please cite the following manuscript if you find TrackSOM useful in your research.

> Givanna H. Putri, Jonathan Chung, Davis N. Edwards, Felix Marsh-Wakefield, Suat Dervish, Irena Koprinska, Nicholas J.C. King, Thomas M. Ashhurst and Mark N. Read. (2022). TrackSOM: mapping immune response dynamics through sequential clustering of time- and disease-course single-cell cytometry data. Cytometry A; doi: https://doi.org/10.1002/cyto.a.24668


The citation in bibtex format:
```
@article{putri2022tracksom,
  title={TrackSOM: Mapping immune response dynamics through clustering of time-course cytometry data},
  author={Putri, Givanna H and Chung, Jonathan and Edwards, Davis N and Marsh-Wakefield, Felix and Koprinska, Irena and Dervish, Suat and King, Nicholas JC and Ashhurst, Thomas M and Read, Mark N},
  journal={Cytometry Part A},
  year={2022},
  publisher={Wiley Online Library}
}
```

R scripts to reproduce evaluations and figures in our paper introducing the `TrackSOM` algorithm are available on [`TrackSOM-evaluations`](https://github.com/ghar1821/TrackSOM-evaluations) repository.

## Support and contribute
We are continuously working to improve TrackSOM and welcome feedbacks.
If you have questions or have encountered issues using TrackSOM, please post an [`issue ticket`](https://github.com/ghar1821/TrackSOM/issues).

## Publications using TrackSOM
TrackSOM was used to analyse COVID-19 temporal data published in Cell Reports Medicine journal:


> Koutsakos, M., Rowntree, L.C., Hensen, L., Chua, B.Y., van de Sandt, C.E., Habel, J.R., Zhang, W., Jia, X., Kedzierski, L., Ashhurst, T.M. and Putri, G.H., 2021. Integrated immune dynamics define correlates of COVID-19 severity and antibody responses. Cell Reports Medicine, p.100208.

R scripts to reproduce evaluations and figures in the paper are available on [`TrackSOM-evaluations`](https://github.com/ghar1821/TrackSOM-evaluations) repository.

## License
TrackSOM is released under the GPL-3.0 lincense, included [here](LICENSE).


## Code of Conduct
  
Please note that the TrackSOM project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.





