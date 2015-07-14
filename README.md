#ncfilter.py

A module that allows to create new netCDF files based on the information in an
already existing one. Has a command-line interface.

The meta-information (everything sans the data proper) is read into an
internal data-structure in the class `NcFilter`. Methods of this class
implement various operations on that meta-data.  New actual data can
be attached to existing or newly created
variables. `NcFilter.write(outfile)` writes the new netCDF-file.

New functionality can be added by deriving from `NcFilter` and
utilizing the helper methods. The comand line interface for new
functionality can be implemented by adding a function to the class
`Commands`

Two commands are implemented exemplarically:

1. The simple `delvar` command deletes a variable.
2. The more complex `compress` packs all suitable variables in a file and saves
   with zip compression.

The implementation of the `compress` command is based on a script by
Dominik Brunner, Empa.

## Usage:

```
usage: ncfilter.py [-h] COMMAND [ARG [ARG ...]] INFILE OUTFILE

Performs operations on a netCDF file.

positional arguments:
  COMMAND     possible commands: ['compress', 'delvar']
  ARG         arguments for commands
                  "compress": ARG = compression_level (int), default=9
                  "delvar": ARG = variable_to_delete (str)
  INFILE      input file
  OUTFILE     output file

optional arguments:
  -h, --help  show this help message and exit

OUTFILE = INFILE will replace the INFILE with the OUTFILE.
```

