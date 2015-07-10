#ncfilter

A class that allows to create new netCDF files based on the information in an
already existing one. Has a command-line interface. Implements packing
and compressing of variables as example

## Status
Work in progress -- don't use.

<!--
# compress_netcdf

Compresses a netcdf file by first converting float32 or float64 type
2D, 3D or 4D fields into short unsigned integers with offset and scaling factor.
Then compress with zlib compression.

## Author:

Dominik Brunner, Empa, Switzerland

## History:

16 June 2015: first implementation

## Usage:

```bash
compress_netcdf.py [-W] -i <infile> -o <outfile>
```
-W : use this option if infile should be replaced by compressed outfile. Option -o <outfile> will be ignored in this case.    
-i <infile>: the orginal non-compressed netcdf file    
-o <outfile>: the compressed netcdf file    
-->
