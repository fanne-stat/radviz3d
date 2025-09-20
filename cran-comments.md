# CRAN Submission Comments

## R CMD check results

0 errors | 0 warnings | 4 notes

* This is a new release.

* Checking package dependencies ... INFO
  Packages suggested but not available for checking: 'MixSim', 'gtools'

  These are suggested packages that provide optional functionality for optimal anchor ordering. The package works correctly without them.

* checking top-level files ... NOTE
  Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed.
  Non-standard files/directories found at top level: 'README.Rmd' 'README_files'

  These are standard development files that are excluded from the built package via .Rbuildignore.

* checking R code for possible problems ... NOTE
  radialvis3d: warning in rgl.spheres(x = 0, y = 0, z = 0, r = radius, color = "grey", add = TRUE, alpha = alpha, depth_mask = FALSE, ...):
  partial argument match of 'r' to 'radius'

  This is a minor warning about partial argument matching in the rgl package function call. The function works correctly.

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.

  This is a system-level issue not related to the package code.

## Test environments
* local macOS install, R 4.5.1
* R CMD check --as-cran

## Description
This package provides functions for creating 3D radial visualizations of multivariate data, extending traditional RadViz techniques to three-dimensional space for enhanced exploration of high-dimensional datasets.