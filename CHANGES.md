# Changes

## v2.0.0

### New features:
* ``zchecker list-found`` to summarize objects already found.
* ``zchecker pccp`` to search the MPC Potential Comet Confirmation Page.

### Enhancements
* Use ``sbsearch`` library.
* Consistent use of ``--stop`` for date ranges.
* ``zchecker clean-eph`` replaced with ``zchecker eph-update --clean``.
* ``zchecker search --vlim`` replaced with the more precise ``--vmax``.
* ``zchecker download-cutouts`` now includes sky reference images.

### Other changes
* All new database format.
* Default magnitude limit is now 25 (V-mag estimates can be _way_ off).
