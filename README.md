## molcat - print molecules to a graphics terminal

molcat is a utility script that simply read a structure file or SMILES and
shows a 2D sketch of the molecule(s) to the terminal. It requires a terminal
that supports the graphics protocol used by, among others, kitty, Ghostty, and
iTerm2.

Example:

![Screenshot of terminal showing efavirenz from SMILES using molcat](https://raw.githubusercontent.com/i-tub/molcat/master/molcat.png)

```
usage: molcat [-h] [-n N] [--all] [--idx | --zidx] [--keeph] [--size-x SIZE_X]
              [--size-y SIZE_Y] [--log-level LOG_LEVEL]
              [file_or_smiles]

Display a 2D sketch of a structure, from a SMILES or a file, to a terminal
that support graphics, such as kitty, Ghostty, and iTerm2.

positional arguments:
  file_or_smiles        structure input file or SMILES strings. If not
                        provided, SMILES will be read from stdin.

options:
  -h, --help            show this help message and exit
  -n N                  index of structure to display. May be a range ('-n
                        1-4')
  --all, -a             show all structures in the file
  --idx, -i             show atom indexes (1-based)
  --zidx, -z            show atom indexes (0-based)
  --keeph, -H           keep all hydrogen atoms
  --size-x SIZE_X, -x SIZE_X
                        X dimension in pixels; default=500
  --size-y SIZE_Y, -y SIZE_Y
                        X dimension in pixels; default is a function of -x
  --log-level LOG_LEVEL
                        RDKit log level; default="FATAL"
```

### Requirements

- Python (tested with 3.11)
- RDKit (tested with 2025.09.6)
