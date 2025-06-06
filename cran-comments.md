# Version 0.1.1

This is the first submission of `ctsframework` to CRAN.

## Comments

Addressed the following comments from CRAN reviewer:

* From: Beni Altmann <benjamin.altmann@wu.ac.at> on 03-Jun-2025:

  - If there are references describing the methods in your package, please add
    these in the description field of your DESCRIPTION file in the form.

    * There are currently no references for this package.

  - Please ensure that your functions do not write by default or in your
    examples/vignettes/tests in the user's home filespace.

    * The code in the vignette was writing to the current working directory.
      This has been corrected. In addition, the vignette has been expanded and
      improved.

  - Please do not set a seed to a specific number within a function.

    * This has been addressed. The RNG seed, which can be specified by the
      user, is stored in an object. When the seed changed, the current RNG
      state is saved and restored afterwards, so there is no impact from the
      user's perspective.

## Test environments

* Local:
  - Windows 11: R 4.5.0 x86_64-w64-mingw32/x64 (64-bit))
* win-builder:
  - R version 4.5.0 (2025-04-11 ucrt) x86_64-w64-mingw32 (64-bit)
    - 1 NOTE:
      * New submission
  - R Under development (unstable) (2025-05-29 r88251 ucrt) x86_64-w64-mingw32 (64-bit)
    - 1 NOTE:
      * New submission

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are currently no downstream dependencies for this package.


