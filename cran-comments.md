## Test environments
* R Under development (unstable) (2022-07-07 r82559)
* Platform: x86_64-pc-linux-gnu (64-bit)
* Running under: Ubuntu 22.04 LTS

## R CMD check results
There were no ERRORs, WARNINGs.

0 errors | 0 warnings | 1 note

* This is a new release.

## Resubmission
This is a resubmission. In this version I have:
* Converted the DESCRIPTION title to title case.
* Converted the DESCRIPTION Authors@R's email address.
* Converted the R/worm_download.R. added mode="wb" to download.file function for 'error reading from connection'.
* Converted the @examples. added if(interactive()){} for 'elapsed time'.

## Resubmission
This is a resubmission. In this version I have:
* Converted the @examples if(interactive()){} to \donttest{}.

## Test environments
* R Under development (unstable) (2022-07-07 r82559)
* Platform: x86_64-pc-linux-gnu (64-bit)
* Running under: Ubuntu 22.04 LTS

## R CMD check results
There were no ERRORs, WARNINGs, Notes.

0 errors | 0 warnings | 0 note

* This is a update.

## Resubmission
* This is a update 0.1.0 -> 0.1.1. 
* issue: Suboptimal check of dist
* https://github.com/rikenbit/WormTensor/issues/24
