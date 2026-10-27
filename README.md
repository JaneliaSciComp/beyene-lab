# Beyene lab analysis

MATLAB analysis code written for the Beyene lab at Janelia, for light-microscopy imaging of dopamine release from dopamine-neuron axons. It was used for Bulumulla et al., "Synapse-like Specializations at Dopamine Release Sites Orchestrate Efficient and Precise Neuromodulatory Signaling", bioRxiv (2024), https://doi.org/10.1101/2024.09.16.613338.

The imaging data are not in this repository. Each script reads its input images from paths set near its top, so change those to point at your own data.

## Channels

- **TDT**: tdTomato, which marks the dopamine-neuron axons
- **ΔF** (or ΔF/F): dopamine released from those axons
- **MAP2**: dendrites
- **BSN**: Bassoon, a presynaptic marker, in the TDT-expressing neurons
- **DRD1**: the D1 dopamine receptor
- **DAPI**: nuclei

## `correlationAnalysis/`: release relative to dendrites and synapses

- `correlations.m` works on a multichannel image (TDT, MAP2 and ΔF):
  - finds release hotspots (local maxima of ΔF within the axons) and compares their distance to MAP2 with that of the axons as a whole, using a two-sample Kolmogorov–Smirnov test
  - measures how hotspots cluster in space, from the pairwise distances between them
  - measures how close the clustered hotspots sit to MAP2
  - computes correlation coefficients between ΔF and the other channels
- `correlationsBetweenBsnDf.m` masks Bassoon and ΔF to the axons and measures how they colocalize: a BSN × ΔF image, a correlation plot and a local correlation image (`getCorrelationImage.m`, `localCorr2.m`).
- `correlationsFromChandima.m` is a version of `correlations.m` from Chandima Bulumulla.
- The remaining files are helpers for heatmaps, confidence intervals and pair indexing.

## `drd1Colocalization/`: receptors relative to axons and nuclei

- `drd1Colocalization.m` thresholds the DAPI, DRD1 and TDT images, takes distance transforms, and bins each channel's signal by distance from the others (`binByDistance.m`, `sumByDistance.m`).
- `drd1ColocalizationRandomShift.m` does the same over 100 trials in which the TDT image is randomly shifted, giving a chance-level baseline.
- In both, `doMask` sets whether DAPI is masked by DRD1. Results are written to a `tifs/` directory.

## `compositeImage/`: synchronized movie

- `compositeMovie.m` renders one frame per time point: the imaging (neurons in red, ΔF/F in a fire lookup table, with a 20 µm scale bar and a time stamp) beside the mean ΔF/F trace with shaded error, drawn up to that time point. The frames are saved as TIFFs; a comment at the end of the script gives the ffmpeg command that joins them into a movie.
- `combineImages.m` applies the color maps, and `plotWithShadedError.m` draws the trace.

## License

BSD 3-Clause; see [LICENSE](LICENSE).
