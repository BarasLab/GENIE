# GENIE-ArtifactFinder

The basic principle is first to simply filter unique variants for which 

- counts in a given panel >= 10

- counts aggregated from all other panels is < 10 and there are >= 3 different panels that cover the variant in question


Of course, we can adjust these thresholds. In particular we need better account for truly different panels from different sites vs different versions of the same "panel" from the same site. And, we will more than likely need some filtering on number of samples per panel.


I also annotated the resulting set with a fisher p-value for the variants that pass this filtering, with the 2 x 2 table of

| | |
| - | - |
| Number of variants called in panel | Number of samples for this panel |
| Number of variants called in other panels with coverage for this variant | Number of samples for panels with coverage for this variant |
