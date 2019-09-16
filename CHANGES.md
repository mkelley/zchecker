# Changes

## v2.4.
* Removed infobits=0 requirement in ZStack and ZPhot
* ZPhot default color correction based on Solontoi et al. 2010.

## v2.4.0
* New `zphot measure` and `zphot plot`.
* New script to plot cutout files with DS9: `ds9-ztf`
* Added ability to download full-frame images for PCCP searches.
* New support for alternate designations.
* Further fix to address duplicate ZTF PID bug.

## v2.3.3
* Smarter use and masking of difference images.
* Fix duplicate ZTF PID bug after ztf-update.
* Enable incremental updates with ztf-update.
* clean-found will remove dates without specifying object.

## v2.3.2
* Download difference images and use for sun-angle image, if available.

## v2.3.1
* Fix missing science image download bug.
* `zstack --clean-missing` to verify and fix local stack file repository.

## v2.3.0
* `zchecker cutout --missing` to verify local cutout file repository.

## v2.2.0
* Cutout sizes can be specified in configuration file and at run time.

## v2.1.1
* Fix download cutouts crash.

## v2.1.0
* Add PCCP cutout downloads.

## v2.0.1
* Catch WCS transformation errors in early commissioning data.

## v2.0.0

### New features:
* ``zchecker list-found`` to summarize objects already found.
* ``zchecker pccp`` to search the MPC Potential Comet Confirmation Page.

### Enhancements
* Use ``sbsearch`` library.
* Optional use of ``--stop`` for date ranges.
* ``zchecker clean-eph`` replaced with ``zchecker eph-update --clean``.
* ``zchecker search --vlim`` replaced with the more precise ``--vmax``.
* ``zchecker download-cutouts`` now includes sky reference images.
* Most scripts can now consider lists of targets.

### Other changes
* All new database format.
* Default magnitude limit is now 25 (V-mag estimates can be _way_ off).
