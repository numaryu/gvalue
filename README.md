# G Value calculation program

Ionization and singlet and triplet excitations.

The original version (S. Sato, Radiation Chemistry 14, 2-19 (1979) [in Japanese]) is written in Fortran 77, 
is rewritten by numaryu into Fortran 2023.

See the [document](docs/gvalue.pdf) for detail.

## References
1. S. Sato, K. Okazaki, S. Ohno, Bull. Chem. Soc. Jpn, 47, 2174 (1974).
2. K. Okazaki, S. Sato, S. Ohno, Bull. Chem. Soc. Jpn, 48, 1411 (1975).

## Build

```
# git clone https://github.com/numaryu/gvalue.git
# cd gvalue
# make
```

Edit make.conf for your environment.

## Run

- Place generated executable (gvalue in src) and input file in workdir. The stem of input file (filename without extension) defines the `runname`.
- Execute the command (`runname` in this case is test).
```
# cd ${workdir}
# ./gvalue test.in 
```
- Results are stored in `runname`.dat and `runname`_degradation.dat.
- Sample GNUPLOT script plot_all.gp (in tools) can be used to generate graphs. Change runname and outfile as you like.
