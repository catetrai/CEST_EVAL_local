CEST evaluation for Bruker, step-by-step
====

----------


Caterina Trainito

https://github.com/catetrai

caterina.trainito@student.uni-tuebingen.de


----------


## Input CEST data ##

The analysis works on 2-D CEST data. If you have acquired more than one slice you will have to run the preprocessing separately for each slice.

For basic Z-spectra analysis, the following scans are necessary and sufficient:

 - **Mz**: one series of saturation images (e.g. 61 images with frequency offsets -7.5 : 0.25 : +7.5 ppm).
 - **M0**: one image acquired without RF saturation pulse.
 
Additional acquisitions that are used for field inhomogeneity corrections and for T1 mapping:

 - **Mz** image series at different B1 strengths: e.g. with B1 = 0.5-3.0 muT. At least two such series are needed for B1-inhomogeneity correction with the 'Z-B1-correction' method (see Windschuh et al., 2015). Note that B1 correction also requires a B1 map obtained with a WASABI sequence.
 - **WASABI**: image series for simultaneous B0- ("WAter Shift") and B1- ("BI") mapping. These maps are used for field inhomogeneity correction (see Schuenke et al., 2016).
 - **T1 mapping sequence**: a series of scans with different
inversion recovery times for T1 mapping. A T1 map is needed for calculating the relaxation-compensated CEST contrast, AREX (see 'Contrasts' below).

## Analysis pipeline ##

For a complete CEST evaluation, the general pipeline is the following: we first process the WASABI series to obtain simultaneous B0- and B1-maps. Then we run B1-correction using multiple B0-corrected Mz series acquired at different B1 values. The corrected Z-spectra are stored in the variable `Z_corrExt`. Finally, we perform pixel-wise Multi-Lorentzian fitting of the corrected Z-spectra. This allows us to evaluate separate CEST peaks (amide, amine, etc.) and calculate CEST contrast images using different metrics.

### Loading data and variables ###
(1) Set the variable `sequenceType` to either 'Mz' (to load one saturation image series and M0), 'WASABI' (to load the WASABI series and M0), or 'T1mapping' (to load Inversion Recovery images for T1 mapping). Note that 'WASABI' uses the same loading routines as 'Mz'.

(2) Choose the respective scan from the study protocol. Indices of the selected protocol entries will be stored in structure `ix`.

(3) Create a `P` structure containing parameters for the fitting functions. Parameters are extracted automatically from Bruker metadata files.

(4) Define an image ROI mask (consisting of 1s and NaNs) called `Segment` that selects pixels for all following analyses. Create it manually with the `make_Segment` UI tool, or load a predefined one.

### Computing Z-spectra ###
The Z-spectrum for each pixel is obtained by normalizing the `Mz_stack` image by `M0_stack`. This is the most basic operation to visualize and quantify CEST effects.
```matlab
Z_uncorr = NORM_ZSTACK(Mz_stack, M0_stack, P, Segment);
```

### Correcting for field inhomogeneities ###

Before computing Z-spectra, ideally you would like to correct the Mz images for inhomogeneities both in the static magnetic field (B0) and in the saturation field (B1). This requires a dB0 map and a relative B1 map, respectively.  If you have acquired a WASABI sequence, you can derive both simultaneously using the WASABI method. Otherwise, you can still perform B0 correction using an internally-computed dB0 map.

#### WASABI evaluation ####
(1) Load WASABI data and do M0-normalization.

(2) Run WASABI fit. This step is implemented by the general fitting function `FIT_3D`, which uses parallel computing (if your Matlab license supports it). A numeric code for the specific WASABI routine is stored in the `P` structure in the field `P.FIT.modelnum`.

(3) From the output of the fit, we calculate the dB0-map (`dB0_stack_ext`) and the relative B1 map (`B1map`). These can be plotted and visualized as images. The order of magnitude of the pixel values will give you a feeling for whether field inhomogeneities are a major concern in your acquisitions.

(4) Make 5-D array `Z_stack` of B0-corrected Z-spectra at different B1 strengths (dimensions: x,y,z,w,B1).

(5) Run Z-B1-correction:
```matlab
Z_stack_corr = Z_B1_correction(Z_stack, B1map, B1_input, B1_output, Segment, 'linear');
```

   `B1_input` is the vector of multiple B1 values that you have acquired. `B1_output` is the B1 value of the Mz image we want to correct (but note that the function can return fitted images for multiple interpolated B1 values).
   
(6) Our final B0- and B1-corrected Z-spectrum is the image at the desired `B1_output` value:
```matlab
Z_corrExt = Z_stack_corr(:, :, :, :, 1);
```

#### B0-correction using internal map ####

(1) Calculate an internal dB0 map. This is given by the x-axis offset of the minimum of the interpolated Z-spectrum from the nominal 0 ppm.
```matlab
dB0_stack_int = MINFIND_SPLINE_3D(Mz_stack, Segment, P);
```
(2)  Then do pixel-wise B0 correction. This simply centers each pixel’s Z-spectrum by shifting it by the calculated x-offset dB0 map):
```matlab
Mz_CORR = B0_CORRECTION(Mz_stack, dB0_stack_int, P, Segment);
```
(3) Finally, compute B0-corrected Z-spectra using `NORM_ZSTACK`.

### T1 mapping ###
Implemented by function `T1eval_levmar`. The fitting is sensitive to starting parameters -- try tweaking them to see if the output improves.

### Multi-Lorentzian fitting of Z-spectra ###

To model CEST effects of interest we fit a sum of Lorentzian line shapes to each pixel’s Z-spectrum. The fit is performed pixel-wise by calling `FIT_3D` with fitting parameters specified in the `P`structure (default is 5-pool model: water, amine, amide, NOE, MT). The function `get_FIT_LABREF` then calculates Zlab and Zref (used for CEST contrasts).

For pH-weighted Ultravist-CEST, use the function:
```matlab
[Zlab, Zref, P, popt] = lorentzianfit_main(Z_corrExt, P, Segment, 'invivo');
```
The ‘invivo’ argument specifies a 6-pool Lorentzian model that includes Ultravist peaks (water, amine, amide, NOE, 5.6ppm, 4.2ppm). For phantom scans, set the argument to ‘ultravist’. This fits a 3-pool model corresponding to Ultravist peaks only (Chen et al., 2014).

`Zlab` is the full n-pool fit (the 'label' image). `Zref` is a structure containing reference images for each pool i (i.e. for each pixel, the sum of all Lorentzians _excluding_ pool i).

### CEST contrasts calculation ###

The CEST effect can be quantified using different metrics:

 - **MTR_asym** _(not implemented)_: Magnetization Transfer Ratio asymmetry. It is simply computed as the difference between each Z-value and the value at the symmetrically opposite spectral location.
 - **MTR_LD**: Magnetization Transfer Ratio, where the asymmetry is calculated as linear difference between Zref and Zlab at each frequency offset. This measure is confounded by water saturation spillover effects (‘dilution’ of the Z-spectrum), especially at higher B1 strengths (Zaiss et al., 2014).
 - **MTR_Rex**: Magnetization Transfer Ratio, where the asymmetry is calculated as the difference of the reciprocal terms. Does spillover- and MT-correction.
 - **AREX** (Apparent Exchange-Dependent Relaxation): MTR_Rex normalized by T1 map. Relaxation-compensated measure.

## Exporting CEST images ##
CEST contrast images can be written to DICOM format. However, DICOM stores data in uint16 type. This means that CEST contrast values (encoded as double-precision floating point, typically smaller than 1) will be automatically re-scaled to unsigned integers. To be able to recover the original values, together with the DICOM file we also save a .mat file ('[filename]_ScalingFactorMap') with the pixel-wise scaling factors. This way, if you load into MATLAB a DICOM CEST image (e.g. after ROI analysis in PMOD) you can convert the pixel values back to the original scale and then do statistics and further analyses on them.

## References ##
Chen, L. Q., Howison, C. M., Jeffery, J. J., Robey, I. F., Kuo, P. H., & Pagel, M. D. (2014). Evaluations of extracellular PH within in vivo tumors using acidocest MRI. Magnetic Resonance in Medicine, 72(5), 1408–1417. https://doi.org/10.1002/mrm.25053

Schuenke, P., Windschuh, J., Roeloffs, V., Ladd, M. E., Bachert, P., & Zaiss, M. (2016). Simultaneous mapping of water shift and B1(WASABI)-Application to field-Inhomogeneity correction of CEST MRI data. Magnetic Resonance in Medicine, 0, 1–10. https://doi.org/10.1002/mrm.26133

Windschuh, J., Zaiss, M., Meissner, J. E., Paech, D., Radbruch, A., Ladd, M. E., & Bachert, P. (2015). Correction of B1-inhomogeneities for relaxation-compensated CEST imaging at 7T. NMR in Biomedicine, 28(5), 529–537. https://doi.org/10.1002/nbm.3283