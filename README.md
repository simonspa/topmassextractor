topmassextractor
================

Extract the top quark mass from event yields or normalized differential cross section. Requires output histograms from the DESY top analysis framework.

# Compilation #

The topmassextractor requires the ROOT libraries and headers to be present. Also, CMake is needed.

To compile the extractor library and the executable, run:

  ```
  $ mkdir build && cd build/
  $ cmake ..
  $ make
  ```


# Usage #

The extraction is invoked by the following command:

```
$ ../bin/extract [options]
```

with the following possible command line arguments:

  * `-t [yield|diffxs|yieldstats|diffxsstats]`: select the type of extraction, either from the total event yield (`yield`) or from the normalized differential cross section (`diffxs`). Also provides the switches for calculating the statistical uncertainty on the systematic variations applied (`yieldstats` and `diffxsstats`).
  * `-v [CRITICAL|ERROR|RESULT|WARNING|INFO|DEBUG|DEBUG2-4]`: select the verbosity level of the extractor. For normal operation, running either `RESULT` or `INFO` should be fine, lower levels may produce a lot of output and slow down the process.
  * `-i [path]`: input path, should be the parent folder of the top mass analysis framework, i.e. the folder containing "preunfolded", "SelectionRoot", "SVD", "UnfoldingResults"
  * `-o [path]`: output path, where all histograms, tables and PDFs will be stored.
  * `-c [ee|emu|mumu|combined]`: select the channel to run on. Of no argument is given, extraction from all channels is performed.
  * `-s`: run on all systematic variations. Without this flag, only the nominal sample is evaluated and the mass extracted.
  * `-d`: do not run a closure test but use real data for the extraction. Only for yields.
  * `--pseudo`: do run on pseudo data for closure (blinded analysis). Only for yields.
  * `-l [filename]`: in addition to screen output, also write logging messages to file given by `filename`.
  * `-m [sample name]`: select the sample to be used as closure data (only has an effect without `-d`). Possible sample names are `MASS_DOWN_{6,3,1}GEV`, `Nominal`, `MASS_UP_{1,3,6}GEV`
  * `-f [token[,token]]`: allows specification of runtime flags. Multiple flags can be given using comma as separator. Do not include any blanks! (I.e. 

# Flags / Flag Tokens #

The following flag tokens for the `-f` command line argument are currently supported:
    
  * `fit | nofit`: Get Chi2 distribution from already fitted bin distributions instead of calculating the Chi2 just at the measurement points. Default is `fit`.
  * `root`: Do create and store histograms and canvases into an output Root file.
  * `pdf`: If flag `root` is set, in addition store all canvases to PDF files into the output directory. This is only active for the nominal extraction, not for systematics.
  * `pdfall`: If flag `root` is set, in addition store all canvases for all systematic variations processed to PDF files into the output directory. This is very slow.
  * `lastbin`: Do only extract from the last bin (most sensitive to the top quark mass) of the histogram instead of the full distribution.
  * `pred | nopred`: Enable disable inclusion of theory prediction uncertainties in the MC statistical error for all extractions. The theory prediction uncertainties taken into account are Q^2 scale and Matching, the errors are calculated by taking the sample difference to nominal and added in quadrature to the statistical error. Default is `nopred`.
  * `norm | nonorm`: Enable/disable normalisation of the total distribution. This flag can be used for event yield and differential cross section. Default is `norm`.
  * `nochlabel`: Do not plot the channel labels

The following flags only apply for the `yield` mode:

  * `bgr | nobgr`: Do not subtract the background. The data is just taken as is, from the MC signal and backgrounds, a "pseudo data" sample is produced including the backgrounds. Default is `bgr`.

The following flags only apply for the `diffxs` mode:

  * `cov | nocov`: enable or disable calculation of bin-to-bin correlations using the covariance matrix from unfolding.
  * `mcstat`: Flag to explicitly exclude/ignore the statistical error on the MC sample in the chi2 calculation. If set, just the data statistical errors (convoluted with whatever has been used for unfolding) are taken into account. This should only be used to evaluate the statistical errors of systematic variations when extracting from differential cross-section, not to extract the mass or any systematic uncertainties.
