# Glacier Response to Subglacial Flowrate Reduction

This repository contains MATLAB scripts for simulating how local water removal from the subglacial drainage system modifies basal drag and glacier speed in an idealized 2D cross section.

The coupled solutions for 2D ice speed and 2D ice temperature are solved following Ortholine v1.0, a coupled thermomechanical, antiplane-strain, free-boundary model developed by Suckale and Elsworth to evaluate shear-margin stability. This repository implements hydrology-informed basal-strength parameterizations to test the response of an idealized mountain glacier to subglacial water removal across different drainage modes. The water is removed either with exact spatial precision, where the offset is y<sub>off</sub> = 0 (on-target flowrate reduction), or it misses a desired target by some non-zero offset, y<sub>off</sub> (off-target flux reduction).

## Files

### Executable scripts

- `v6Neel_2Dcvx_dim_OnTarget.m`  
  Runs the on-target flowrate-reduction experiments. The intervention is colocated with the drainage feature whose effective pressure and basal strength are being modified.

- `v6Neel_2Dcvx_dim_2SideCanals_OffTarget.m`  
  Runs the off-target flux-reduction experiments. The script imposes a localized basal-strength perturbation at a prescribed off-target suction location in the water-film region, while retaining a central canal and two side canals.

### Input/output folders

- `InputParameterFiles/`  
  Contains `.mat` files for the on-target drainage-mode cases.

- `ResultsMatFiles/`  
  Stores saved output `.mat` files.

- `lib/`  
  Contains mesh-generation/helper functions, including DistMesh-style utilities.

- `cbrewer2/`  
  Contains the ColorBrewer plotting utility.

## Dependencies

Required:

1. MATLAB
2. [CVX for MATLAB](http://cvxr.com/cvx/)
3. DistMesh-style mesh utilities, including `distmesh2d`, `dpoly`, `huniform`, and `boundedges`
4. `cbrewer2`
5. `iceColorMap.mat` for temperature plots

The scripts currently call:

```matlab
gpuD = gpuDevice;
```

If no GPU is available, comment out this line. The rest of the scripts do not explicitly use GPU arrays.

## Usage

### 1. Set up folders

Before running, make sure the repository contains:

```text
InputParameterFiles/
ResultsMatFiles/
lib/
cbrewer2/
iceColorMap.mat
```

In MATLAB, make sure CVX is installed and initialized with:

```matlab
cvx_setup
```

### 2. Run an on-target case

Open:

```text
v6Neel_2Dcvx_dim_OnTarget.m
```

Choose one drainage-mode case by uncommenting exactly one `load(...)` command in the input-file section. For example:

```matlab
load("InputParameterFiles\Canal_50kPa_Baseline.mat");
```

Then run:

```matlab
v6Neel_2Dcvx_dim_OnTarget
```

Available on-target cases include sedimentary canals, R-channels, linked cavities, and water films. Each case is controlled by its corresponding `.mat` file in `InputParameterFiles/`.

Before each run, update the output filename at the bottom of the script so it matches the selected case:

```matlab
filepath = "ResultsMatFiles\";
filename = filepath+"HighRes_"+"LinkedCavity_"+"500"+"kPa.mat";
save(filename)
```

### 3. Run an off-target case

Open:

```text
v6Neel_2Dcvx_dim_2SideCanals_OffTarget.m
```

Set the suction location near the top of the script:

```matlab
Xsuction = 1100; % [m]
```

Suggested values currently listed in the script are:

```text
1010 m, 1020 m, 1050 m, 1100 m
```

Then run:

```matlab
v6Neel_2Dcvx_dim_2SideCanals_OffTarget
```

The off-target script currently uses a canal baseline with `N_canal = 50 kPa`, a localized basal-strength perturbation of `N_wfblip = 1333 Pa`, and a Gaussian length scale `delta = 10 m`. The output file is saved automatically using the selected `Xsuction` value:

```matlab
filename = filepath+"OffTarget_Canal50kPa_Blip1333Pa_XSuction"+string(Xsuction)+".mat";
```

If you change `Xsuction`, also update any manually named diagnostic variables in the AUC/basal-drag saving section if you want the variable names to match the case.

## Output

Each script generates plots of:

1. 2D ice speed field
2. Ice surface speed profile
3. 2D temperature field
4. Basal drag profile

The scripts also save the final MATLAB workspace as a `.mat` file in `ResultsMatFiles/`.

## Notes

- The on-target script reads hydrology-informed basal-strength functions from `.mat` files.
- The off-target script defines the mixed-bed hydrology and off-target perturbation directly inside the script.
- Windows-style paths are currently used. For better cross-platform compatibility, replace paths such as

```matlab
load("InputParameterFiles\Canal_50kPa_Baseline.mat");
```

with

```matlab
load(fullfile("InputParameterFiles", "Canal_50kPa_Baseline.mat"));
```

## Citation and acknowledgment

The coupled thermomechanical/free-boundary model is adapted from Ortholine v1.0:

```text
Suckale, J. and Elsworth, C. W.: An antiplane strain model for evaluating shear-margin stability (Ortholine v1.0), EGUsphere [preprint], https://doi.org/10.5194/egusphere-2026-67, 2026.
```

Code repository for the original model:

```text
https://github.com/coopere/InstituteIceStream2D
```

For use of this modified hydrology/intervention code, please also cite the associated basal-anchoring manuscript once available.

## License

Add a license before public release. Common options for academic research code include MIT, BSD-3-Clause, GPL-3.0, or an institutional license.

## Contact

For questions, please contact the repository maintainer.
