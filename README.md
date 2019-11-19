# DEM2GRD

This project is no longer maintainted.

Please see the new project @ https://github.com/mattbilskie/PyDEM2GRD.

This program assigns ADCIRC mesh nodes an elevation based on a digital elevation model (DEM).

## Citation

M.V. Bilskie, S.C. Hagen (2013) Topographic Accuracy Assessment of Bare Earth lidar-derived Unstructured Meshes. Advances in Water Resources, 52, 165-177, doi:10.1016/j.advwatres.2012.09.003.

## Getting Started

### Prerequisites

cmake version 2.8.12 or higher

### Installing

```
mkdir build; cd build
ccmake ..
make
make install
```

## Running the Program

### Inputs

This program requires three inputs. The first is a general control input file that contains required information/parameters to run the program (input.inp). The other inputs are an ADCIRC mesh and DEM (in flt/hdr format).

#### Control Input File

The contents are as follows:

* Line 1: Header Line
* Line 2: Header Line
* Line 3: Name of flagged ADCIRC mesh file to be interpolated
* Line 4: Coordinate system (0 for Cartesian and 1 for Geographic)
* Line 5: Z-multiplication factor
* Line 6: Output ADCIRC mesh wither interpolated values
* Line 7: Number of flt/hdr rasters
* Line 8+: Raster file names (without the file extension)

Notes on the flagged ADCIRC mesh:

Only flagged nodes (less than or equal to -1000) will be included for interpolation. A flag value of -1001 uses the method in Bilskie & Hagen. A flag value of -1002 will use the methods in Bilskie & Hagen and multiply the result by 2, a value of -1003 will multiple by 3, and so on. This essentially increases the control volume and is a method to smooth the elevations. Including a “1” in the hundredths place will only interpolate “wet” value within the control volume (i.e. -1101). This is useful for interpolating within rivers and channels when you do not want to include DEM cells on the land in the averaging method.

## Disclaimer

This program is to be used as-is, and the developers do not provide warranty of any kind. The software may not be error free and may not be appropriate for all projects or decisions. The developers have undergone basic checks, but has not undergone a formal Quality Assurance of any kind. Please report any bugs to mbilsk3@lsu.edu
