## Test environments
* local OS X install, R 3.2.2
* ubuntu 12.04 (on travis-ci), R 3.2.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 3 notes

**checking CRAN incoming feasibility ... NOTE
Maintainer: 'Thomas Clarke <tclarke@jcvi.org>'

New submission

** running examples for arch 'i386' ... [407s] NOTE
Examples with CPU or elapsed time > 10s
                   user system elapsed
ggrasp.cluster   102.60  13.31  115.94
ggrasp.recluster  87.94  10.53   98.47
print.ggrasp      86.36   9.97   96.36
ggrasp.write      82.95   9.87   92.82
** running examples for arch 'x64' ... [399s] NOTE
Examples with CPU or elapsed time > 10s
                   user system elapsed
ggrasp.cluster   105.28  14.58  119.85
ggrasp.recluster  80.09  11.77   91.87
print.ggrasp      79.68  12.12   91.81
ggrasp.write      80.72  10.80   91.51

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

* FAILURE SUMMARY

* All revdep maintainers were notified of the release on RELEASE DATE.
