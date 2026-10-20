# JEWEL 2.6.0 → C++23 porting rules

Goal: a **faithful, line-by-line translation** of the JEWEL 2.6.0 Fortran core
(`jewel-2.6.0.f`, `medium-simple.f`, `medium-vac.f`, `meix.f`) into C++23 that
produces **bit-identical events** to the gfortran-compiled original for the same
random-number seed (both built with `-O2 -ffp-contract=off`). The modified
PYTHIA 6.4.25 stays Fortran and is linked in; JEWEL talks to it via the
common blocks in `include/jewel/pythia6.hpp`.

This is a translation, not a rewrite. Do **not** "improve" physics, change
algorithms, reorder floating-point operations, fix apparent bugs, or skip code
paths. If you find a genuine bug in the original, keep the original behaviour
and report it (comment `// NOTE(port): ...`).

Reference sources: `../src_orig/jewel-2.6.0/*.f`.

## Layout

| header | content |
|---|---|
| `jewel/farray.hpp` | `FArray` – Fortran-layout arrays, 1-based, `a(i,j)` |
| `jewel/pythia6.hpp` | PYTHIA commons (`jewel::py::pyjets::K`, …) and wrappers `py::pyr()`, `py::pyp(i,j)`, `py::pyrobo(...)`, … |
| `jewel/commons.hpp` | every JEWEL COMMON block as `namespace jewel::cb::<block>` |
| `jewel/routines.hpp` | prototypes of all core routines (do not change signatures on your own, see below) |
| `jewel/medium.hpp` | `Medium` interface; core calls `medium().gettemp(x,y,z,t)` etc. |
| `jewel/fortran_compat.hpp` | `r4()`, `ipow<N>()`, `fsign`, `nint`, `fint` |
| `jewel/io.hpp` | `io::log(...)` (JEWEL log = Fortran unit), `io::hepmc()` (FILE*) |

All code lives in `namespace jewel`. File header comment: `// JEWEL 2.6.0 C++ port -- <what>. Translated from <file>:<lines>.`
Keep the original `C--` comments (as `//` comments) – they are the only documentation.

## Common blocks

A Fortran routine that declares `COMMON/PARAM/...` and `COMMON/PYJETS/...` becomes

```cpp
void makebranch(const int& l) {
  using namespace cb::param;
  using namespace py::pyjets;
  ...
```

Names: JEWEL common members are **lower case** (`nscatcen`, `mv(i,j)`, `trip(i)`);
PYTHIA members are **upper case** (`N`, `K(i,j)`, `P(i,j)`, `V(i,j)`, `MSTU(i)`, `PARJ(i)`, `MSTP`, `PARP`, `MSTI`, `PARI`, `CKIN`, `MSUB`, `MSEL`, `KFIN(i,j)`, `MDCY`, `MDME`, `MRPY`).
Fortran is case-insensitive: normalise every reference accordingly.
Local variables: lower case, declared at first use where that is unambiguous, otherwise at the top.

## Arguments & calls

* Signatures are in `routines.hpp`: `T&` for arguments the routine writes, `const T&` otherwise.
  This preserves Fortran pass-by-reference (including aliasing with common variables).
* If you find that an argument marked `const` is actually written (directly or by a callee), or an
  output is marked wrongly, **do not edit routines.hpp** – use the correct signature in your .cpp, and list the
  change in your final report; the integrator merges header changes.
* Calling with a literal/expression for a `T&` parameter: make a local temporary, exactly as Fortran does.
* `CHARACTER` values → `std::string` (compare with `==`; Fortran blank padding is irrelevant for the values
  used here: 'Q','G','QQ','GC','p+','n0','PPJJ', ...). `GETSSCAT` takes `CHARACTER*1` → `char`; callers pass `s[0]`.
* `LOGICAL` → `bool`.

## Numeric rules (bit-exactness!)

1. **Single-precision literals.** A Fortran literal without `d` exponent (`0.4`, `1.e-5`, `0.3333333`) is REAL*4.
   If it is not exactly representable as float→double (0.5, 2., 1.5, 10. are; 0.1, 0.4, 0.3333, 1.e-5 are not)
   write `r4(0.4f)`. When in doubt use `r4(...f)` – it is always correct. Integers-valued literals (`2.`, `100.`) may be
   written as `2.0`.
2. **Pure single-precision sub-expressions** are evaluated in single precision: `2./3.` → `r4(2.f/3.f)`,
   `0.1*0.2` → `r4(0.1f*0.2f)`. Only when *both* operands are REAL*4 literals/variables.
3. Double-precision literals `1.d-4`, `3.141592653589793d0` → plain C++ double literals.
4. **Integer powers:** `x**2` → `ipow<2>(x)` (or `x*x` for a plain variable), `x**4` → `ipow<4>(x)`, `x**6` → `ipow<6>(x)`.
   Never `std::pow` for integer exponents. Integer base & exponent: `ipowi<2>(n)`.
   Real exponent: `x**0.25` → `std::pow(x, 0.25)`, `x**0.3333333` → `std::pow(x, r4(0.3333333f))`,
   `x**(1./3.)` → `std::pow(x, r4(1.f/3.f))`, `x**2.` → `x*x`.
5. **Mixed integer arithmetic:** `1/2` is integer division (0). `n*1.d0/m` evaluates left to right.
   `DBLE(n)`, `1.d0*n` → `static_cast<double>(n)` / `1.0*n`.
6. Keep the **evaluation order and parenthesisation exactly** as written; `a*b*c` is `(a*b)*c` in both
   languages. Do not factor, hoist or simplify arithmetic.
7. Intrinsics: `SQRT, EXP, LOG, LOG10, SIN, COS, TAN, ATAN, ATAN2, ACOS, ASIN, COSH, SINH, TANH, ABS` → `std::` counterparts;
   `SIGN(a,b)` → `fsign(a,b)`; `NINT` → `nint`; `INT` → `fint` (truncation); `MOD(i,j)` (ints) → `i % j`;
   `MAX/MIN` (≥2 args) → nested `std::max/std::min` preserving argument order; `DBLE/DFLOAT` → cast.
7a. **SIN/COS pairs.** gfortran merges `SIN(X)` and `COS(X)` of the same value inside a routine into a single
   `cexp(0+iX)` call; clang would instead merge them into `__sincos_stret`. Both can differ from separate
   `sin`/`cos` in the last bit. The core is compiled with `-fno-builtin-sin -fno-builtin-cos`, and exactly the
   sites that gfortran merged (determined from the reference objects: MAKESPLITTING, MAKEINSPLIT, GETQVEC,
   DOKINEMATICS ×3, the first split block of DOCOHSCATTERING, GETSCATTERER of the medium) use
   `gf_sin`/`gf_cos` (fortran_compat.hpp). Check with
   `objdump -d -r <obj> | grep -E "RELOC_BRANCH26.*_(cexp|sin|cos)$"` after changes.
8. `.AND.`/`.OR.` → `&&`/`||` (gfortran -O2 short-circuits left to right, same as C++). Keep operand order –
   this matters because operands can call `pyr()`.
9. **Random numbers:** every `PYR(0)` → `py::pyr(0)`. The *number and order* of PYR calls must be identical to
   the Fortran; never cache, merge, hoist or reorder them.

## Control flow

* **DO loops:** the trip count is fixed on entry; the loop variable is not affected by later changes of the bound
  (JEWEL often adds lines to the event record, changing `N`, inside `DO I=1,N`). Translate
  `DO 10 I=A,B` as
  `for (int i = a, i_end = b; i <= i_end; ++i)`. If the loop index is used after the loop, declare it
  outside; after normal completion it equals `i_end+1` (for step 1). Negative/other steps: compute the trip count
  `max((b-a+s)/s, 0)` like Fortran. Never let a loop body modify the loop index unless the Fortran does.
* `GOTO`: prefer structured translation (`while`, `break`, `continue`, early `return`) when it is obviously
  equivalent; otherwise use C++ `goto` with labels named `l<number>:`. Correctness over elegance. Note that a
  C++ `goto` cannot jump over initialised declarations – declare those variables before the first label.
* `RETURN` → `return` (for functions return the current value of the result variable; keep a local result
  variable named like the function when the Fortran assigns it in several places).
* `CALL EXIT(1)` / `STOP` → `std::exit(1)` / `std::exit(0)`.

## Locals, SAVE, DATA

* `DATA x/…/` makes a local **static** (SAVE). If the routine never modifies it, use `constexpr`/`const`;
  if it may modify it, `static`.
* Local variables without SAVE/DATA are, in gfortran, uninitialised stack variables — except arrays larger than
  64 kB, which are static. If a local may be read before it is assigned (on some path), mark it
  `// NOTE(port): possibly uninitialised in Fortran` and initialise it to 0; report it.
* Local arrays → `FArray` (e.g. `F1<double,4> q;`) to keep 1-based indexing.

## I/O

* `write(logfid,*) 'text', x` → `io::log(" text {}", x)`; numbers may be formatted with `{}`
  (exact layout of log lines is not important, content is). `write(*,*)` → `std::println(...)`.
* `call flush(logfid)` / `CALL FLUSH(3)` → `io::log_flush()`.
* HepMC output: `std::fprintf(io::hepmc(), ...)`; formats must reproduce the Fortran HepMC IO_GenEvent text
  (numbers may be printed with enough digits to round-trip; the file is compared numerically).
* Table files (`splitint.dat`, `pdfs.dat`, `xsecs.dat`): C++ streams; read whitespace-separated numbers (Fortran
  list-directed files must be readable), write with `std::format("{:.17e}")` so Fortran can read them back.

## Style (C++23)

`clang++ -std=c++23 -Wall -Wextra`. Use `std::string`, `std::format`/`std::println`, `constexpr`, `[[maybe_unused]]`,
range-for where the Fortran loop is a plain iteration with no index use after the loop. No macros, no `using
namespace std` at file scope. Keep one C++ function per Fortran routine, in source order, with the Fortran name.
Do not introduce classes for the core – the state is global by design (it is shared with PYTHIA).

## Checking your work

Your files must compile standalone: `clang++ -std=c++23 -fsyntax-only -Wall -Wextra -Iinclude src/<file>.cpp`.
Re-read your translation side by side with the Fortran once more at the end, specifically for: missed
single-precision literals, integer division, loop bounds that change inside loops, PYR call order,
`.EQ.` comparisons on strings, and GOTO targets.
