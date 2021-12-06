# The Commander parameter file

## Format

Commander is controlled through a standard parameter file in which each line contains a parameter on the form
```
PARAMETER_NAME = value
```
In case there are multiple lines with the same PARAMETER_NAME, only the first occurance is used. Comments are indicated with `#`, and blank lines are allowed.

Parameter files may be nested through an include statement, `@INCLUDE parfile2.txt`. This will effectively insert all parameters in parfile2.txt at the current position in the original parameter file. Parameters may also be specified through the command line interface as `--PARAMETER_NAME=value`.

## Parameter types

In general, most Commander parameters may be classified in terms of three groups:
1. Infrastructure parameters; these control the behaviour of the algorithms and IO
2. Data set parameters; these define the properties of each data set
3. Model parameters; these define the properties of each sky component

Specification of all currently supported parameters is provided in the next section. Note, however, that Commander is in constant development, and sometimes new features are added before the documentation has been updated. If you find that the code crashes while complaining about a missing parameter that is not documented, please notify both developers and users at `forums.beyondplanck.science`.

## Error checking and debugging

The entire parameter file is parsed at the beginning of the run, and stored in-memory in a Fortran type called `cpar` ("Commander parameters"). Active filenames are validated with respect to existence, but not content. Many types of parameter file errors are therefore automatically detected immediately, but not all. If you find that the code crashes with something that looks like a parameter error, typical things to check are the following:
- Is a given parameter of the correct expected type? (Hint: Check which line causes the crash, and look it up in `comm_param_mod.f90`; that will often reveal which parameter causes the crash.)
- Does a given file contain the expected data type?
- If it is an ASCII input file, does it have the correct format?

