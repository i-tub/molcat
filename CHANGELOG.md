# Changelog

## [0.3.0] - 2026-04-30

### Improvements

- Determine window size automatically.
- Optionally generate a single image with all of the structures in multiple
  columns.

### Fixes

- Use chunking under kitty. It's necessary when the images get large enough.

## [0.2.0] - 2026-04-26

### Improvements

- Copy to clipboard using kitty protocol.
- Write PNG to file.
- Show only structures with a minimum number of atoms.
- Show stereo labels (R/S/E/Z).
- Recognize quoted SMILES when parsing plain text.

### Fixes

- Loss of cis/trans stereochemistry.

## [0.1.0] - 2026-04-21

Initial release.
