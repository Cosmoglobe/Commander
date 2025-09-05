This is the directory for development of the low-level re-analysis of FIRAS data. It contains the following:

- `notes`: Information about the FIRAS instrument and data.
- `original_pipeline`: The files used to generate the original FIRAS data products.
- `reference`: Reference datasets (e.g. sampling rate, etc).
- `src`: Source code for the FIRAS re-analysis. To run, use `python -m {DIR}.{MOD}`, where `{DIR}` and `{MOD}` correspond to what you want to run, e.g. use `python -m pipeline.main` to run the pipeline.

Products of the current analysis can be found here.