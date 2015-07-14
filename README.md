#ncfilter

A class that allows to create new netCDF files based on the information in an
already existing one. Has a command-line interface. Implements packing
and compressing of variables as example (`compress`) and deletion of a variable (`delvar`).

The implementation of the `compress` command is based on a script by Dominik Brunner, Empa. 



## Usage:

```bash
usage: ncfilter.py [-h] COMMAND [ARG [ARG ...]] INFILE OUTFILE

Performs operations on a netCDF file.

positional arguments:
  COMMAND     possible commands: ['compress', 'delvar']
  ARG         arguments for command
                  compress compression_level (int)
                  delvar variable_to_delete (str)
  INFILE      input file
  OUTFILE     output file

optional arguments:
  -h, --help  show this help message and exit

OUTFILE = INFILE will replace the INFILE with the OUTFILE.
```

