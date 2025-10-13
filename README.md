# Zeolite-analyser code

This code measure some geometric features for an overall database. It enables the assessment of a multitude of geometric characteristics of the frameworks, and it reads CIF files of the periodic structures with P1 space group and SiO2 compositions.  It is a new version of the https://github.com/salrodgom/GeoZeoForCIFFileDatabase code. Some new functionalities have been added.

Compilation
===========
cd src<br>
make 

to clean:<br>
make clean

Citation
========
SRG. Balestra, N. Rodríguez-Sánchez, D.Mena-Torres, AR. Ruiz-Salvador, Structural Features and Zeolite Stability: A Linearized Equation Approach, Cryst. Growth Des. , Feb., 2024, ![10.1021/acs.cgd.3c00893](https://zenodo.org/badge/DOI/10.1021/acs.cgd.3c00893.svg)

Output files
============
For each processed framework the executable emits a collection of plain-text files whose names follow the pattern `<structure>_<property>.dat`. They provide the raw geometric descriptors that the program later summarises in console histograms and in follow-up analysis modules. The table below gives a quick reference for new users.

| File pattern | What it contains | Typical use |
|--------------|------------------|-------------|
| `<structure>_sio.dat` | Every Si–O bond length detected in the structure (one per line). | Input for Si–O bond-length statistics and histograms. |
| `<structure>_oo.dat` | O–O distances between oxygen atoms sharing the same silicon in a tetrahedron. | Deriving averages and dispersion of O–O link lengths. |
| `<structure>_osio.dat` | Internal O–Si–O angles within each SiO₄ tetrahedron, in degrees. | Assessing tetrahedral shape and distortion. |
| `<structure>_siosi.dat` | Si–O–Si angles that describe how tetrahedra connect through bridging oxygens. | Characterising connectivity between neighbouring tetrahedra. |
| `<structure>_ooo.dat` | O–O–O angles located inside channels identified during the analysis. | Inspecting channel geometry. |
| `<structure>_sisisi.dat` | Consecutive Si–Si–Si angles along chains of tetrahedra. | Measuring opening angles along silicon chains. |
| `<structure>_sisi.dat` | Si–Si distances between silicon atoms linked through a single oxygen. | Calculating Si–Si bond statistics, NMR parameters, and staggering angles. |
| `<structure>_Bonds.dat` | Detailed log of every atom triplet considered for angle (“bend”) calculations, including indices, labels, and the resulting angle. | Full traceability of the angular data set. |
| `<structure>_Q.dat` | Tetrahedral order parameter \(q\) for each silicon atom, following the Chau & Hardwick metric. | Quantifying the degree of tetrahedral distortion per silicon. |
| `<structure>_Staggering.dat` | Minimum torsion angles (“first- and second-order staggering”) from O–Si–Si–O dihedrals for each connected Si pair. | Evaluating torsional staggering between neighbouring tetrahedra. |
| `<structure>_29SiNMR.dat` | Estimated ^29Si chemical shifts obtained from four literature correlations (in-house, Fyfe, Dawson, Thomas) plus the silicon index. | Predicting solid-state ^29Si NMR signatures. |
| `<structure>_si.dat` | Summary per silicon atom with its \(q\) value and statistical moments (mean, deviation, skewness, kurtosis) of surrounding Si–O, Si–Si, O–Si–O, Si–O–Si, and O–O data. | Compact descriptor set for each tetrahedral site. |
| Standard output (`Histo: …`) | Console report produced after reopening the `.dat` files to build histograms, including mean, min/max, deviation, skewness, and kurtosis. | Quick on-screen diagnostics; the same descriptors feed later thermodynamic correlations. |
