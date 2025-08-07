# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)

## [Unreleased]
### Fixed
- Updated to wf-template v5.6.2, fixing:
    - Sequence summary read length N50 incorrectly displayed minimum read length, it now correctly shows the N50.
    - Sequence summary component alignment and coverage plots failed to plot under some conditions.
- `store_dir` parameter format incorrectly declared in the schema. This does not affect this workflow as it does not use the storeDir directive and has been changed to maintain compliance with our latest testing standard.
- partial ssAAV counts missing par_icg subtype counts.
### Changed
- Updated to wf-template v5.6.2, changing:
    - Reduce verbosity of debug logging from fastcat which can occasionally occlude errors found in FASTQ files during ingress.
    - Log banner art to say "EPI2ME" instead of "EPI2ME Labs" to match current branding. This has no effect on the workflow outputs.
    - pre-commit configuration to resolve an internal dependency problem with flake8. This has no effect on the workflow.
    - Unexpected workflow parameters now cause the workflow to fail.
    - Sample sheets must not contain an alias that starts with "barcode".
- Updated project description.
- Minor decrease to some memory directives to avoid “Process requirement exceeds available memory” errors when running in WSL.

## [v1.2.0]
### Updated
- Diagrams in the README describing AAV genome types. 
### Fixed
- Rounding issues in genome structures table.
- ITR locations can be specified in a BED file with `--transgene_bed`.
### Added
- Ability to load any number of non-transgene plasmid reference files in a directory with `--non_transgene_refs`.

## [v1.1.2]
### Updated 
- Updated Ezcharts to v0.11.2.

## [v1.1.1]
## Added
- Publish report missing from v1.0.3.
### Fixed
- Automated basecaller detection not finding a basecaller model.

## [v1.1.0]
### Changed
- Increase memory allocated to the medaka process.
- Updated Medaka to v1.12.0.
### Removed
- Options for specifying Medaka model (`--medaka_model` and `--basecaller_cfg`) as this is now inferred automatically from the input data.
### Added
- Tagging of BAMs with `AV:Z` specifying AAV genome type.
- Option `--output_genometype_bams` to output BAMs by assigned genome type.
- Support for creating IGV JSON config files.
- Automatic detection of basecalling model and medaka model selection.
- `--override_basecaller_cfg` parameter for cases where automatic basecall model detection fails or users wish to override the automatic choice.

## [v1.0.3]
### Fixed
- Datatype inference error during CSV loading. 

## [v1.0.2]
### Fixed 
- Missing data when using reference files containing FASTA header descriptions.

## [v1.0.1]
### Fixed
- `<img>` tags in the docs.

## [v1.0.0]
### Added
- Added memory and CPU directives to all processes.

### Changed
- Updated docs.

## [v0.0.2]
### Fixed
- Related protocols link

## [v0.0.1]
### Added
- First release of wf-aav-qc
