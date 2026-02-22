# edabrt: Electrostatic Deflector Aberrations

Computes first and second order aberrations of an electrostatic deflector
in the horizontal x-a plane using exact analytic formulas.

## Build

```
gcc -Wall -O2 -o edabrt edabrt.c -lm
```

Or using the Makefile:

```
make
```

## Usage

**Interactive mode** -- run without arguments and follow the prompts:

```
./edabrt
```

**Command-line mode** -- supply all four deflector parameters:

```
./edabrt <r> <ang> <n1> <n2>
```

- `r` -- reference orbit radius (m)
- `ang` -- central angle spanning the deflector (degrees)
- `n1` -- first order electrostatic field inhomogeneity coefficient
- `n2` -- second order electrostatic field inhomogeneity coefficient

## Reference

E. Valetov and M. Berz, *Derivation of Analytic Formulas for Electrostatic
Deflector Aberrations, and Comparison with the Code COSY INFINITY*,
MSUHEP-180212, Michigan State University (2018).

## License

MIT -- see [LICENSE.md](LICENSE.md).
