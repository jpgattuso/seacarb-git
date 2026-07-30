# v3.4.1
## Submission summary

seacarb 3.4.1 is a patch release. It corrects the pH scale of the
Perez and Fraga (1987) hydrogen fluoride dissociation constant (Kf,
argument kf = "pf"), which was previously treated as total-scale and
is in fact free-scale (confirmed by F. Perez, 2026). The unwarranted
total-to-free conversion has been removed, the unused Ks_p0 argument
has been dropped from Kf(), and the internal callers of Kf() and the
Kf documentation have been updated accordingly. The Dickson and Riley
(1979) option (kf = "dg") is unaffected.

## Test environments

* local: Ubuntu 24.04, R 4.4.1
* GitHub Actions: ubuntu-latest (release, devel), macos-latest
  (release), windows-latest (release) — all passing

## R CMD check results

Status: 2 NOTEs, 0 errors, 0 warnings.

* checking for future file timestamps ... NOTE
  unable to verify current time

  This is a network note: the check could not reach a time server to
  verify timestamps. It is unrelated to the package.

* checking HTML version of manual ... NOTE
  Skipping checking math rendering: package 'V8' unavailable

  V8 is not installed in the local check environment. This note
  concerns the check environment, not the package, and does not appear
  on systems where V8 is available.

## Reverse dependencies

seacarb has three reverse dependencies on CRAN (ecolMod, marelac,
respirometry). We checked all three against this version and found no
problems. The kf = "pf" default now returns slightly different Kf
values (altering downstream carbonate-system results only in the
fourth significant figure), and none of the three relies on the
previous values.

