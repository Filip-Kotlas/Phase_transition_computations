# Thermodynamic databases

CALPHAD databases in the TDB format, read by the solver (`include/io/TdbDatabase.hpp`). A model's config names one in `parameters.database` (path relative to the repo root), together with the database's names of the phases it uses, for example:

```json
"database": "databases/NbZr-91Fer-LB.tdb", "alpha_phase": "HCP_A3", "beta_phase": "BCC_A2"
```

Every run copies the database it read into its `info/setup/`.

## Source

All files come from <https://phasediagrams.org/>. Each has a page at `https://phasediagrams.org/phase-diagram/<file name>`. The site plots the phase diagrams with pycalphad. Most of its datasets, including all of these, come from Bengt Hallstedt's collection:

> B. Hallstedt, "The SGTE collection of binary datasets", *Calphad* **89** (2025) 102833, doi:[10.1016/j.calphad.2025.102833](https://doi.org/10.1016/j.calphad.2025.102833).

NbZr-19Laf was downloaded on 2026-11-13 (again, unchanged, on 2026-11-25). The others were downloaded on 2026-11-25. The files are kept as downloaded. A database that has to be changed gets a new file name, because runs record only the name in `parameters.txt` (the copy in `info/setup/` keeps the content).

**Terms of use:** phasediagrams.org states no licence and no terms of use for the files (checked 2026-11-25: home, About, Upload and the files' pages). The files are kept here with their origin and references. Ask the site's maintainer or B. Hallstedt before using them beyond this research code.

## Files

When a database is used in results, cite the original assessment and the collection above. Each file gives its references in its header and in its `LIST-OF-REFERENCE`. The unary data of the elements, which every file contains, come from A. T. Dinsdale, *Calphad* **15** (1991) 317–425, and, for H₂ and the gases, from the SGTE substance database.

| File | System | Assessment | Dataset by Hallstedt |
|---|---|---|---|
| `NbZr-91Fer-LB.tdb` | Nb–Zr | A. Fernández Guillermet, *Z. Metallkd.* **82** (1991) 478–487, doi:[10.1515/ijmr-1991-820609](https://doi.org/10.1515/ijmr-1991-820609) (also in Landolt–Börnstein Vol. 4) | 2011-04-02, modified 2020-10-13 |
| `NbZr-19Laf.tdb` | Nb–Zr | P. Lafaye, C. Toffolon-Masclet, J.-C. Crivello, J.-M. Joubert, *Calphad* **64** (2019) 43–54 | 2020-10-13 |
| `ZrH-97Ans-LB.tdb` | Zr–H | I. Ansara, N. Dupin, unpublished (1997), included in Landolt–Börnstein Vol. 3 | 2012-10-14, modified 2021-02-07 and 2024-08-03 |
| `ZrH-99Dup.tdb` | Zr–H | N. Dupin, I. Ansara, C. Servant, C. Toffolon, C. Lemaignan, J. C. Brachet, *J. Nucl. Mater.* **275** (1999) 287–295, doi:[10.1016/S0022-3115(99)00125-7](https://doi.org/10.1016/S0022-3115(99)00125-7) | 2024-08-04 |
| `ZrH-00Koe.tdb` | Zr–H | E. Königsberger, G. Eriksson, W. A. Oates, *J. Alloys Compd.* **299** (2000) 148–152 | 2024-08-03 |
| `ZrH-07Uki.tdb` | Zr–H | S. Ukita, H. Ohtani, M. Hasebe, *J. Jpn. Inst. Met.* **71** (2007) 721–729 | 2024-08-04 |
| `ZrH-12Zho.tdb` | Zr–H | Y. Zhong, D. D. MacDonald, *J. Nucl. Mater.* **423** (2012) 87–92 | 2024-08-20 |

## Notes on the data

These notes come from the files' own comments and from checks of the reader against pycalphad.

- **NbZr-91Fer-LB** has the same coefficients as the hand-written Zr–Nb thermodynamics used before the TDB reader (commit 458ed20). The current Zr–Nb configs use it.
- **ZrH-97Ans-LB** and **ZrH-99Dup** hold the same assessment. In 97Ans-LB, however, the δ hydride is the separate phase `ZRH2_C1`, and its `FCC_A1` has no Zr:H end member. The two files also differ in `G(ZRH2_L2B,ZR:VA)`, which Hallstedt calls "probably a mistake in the 97Ans-LB dataset".
- **ZrH-12Zho** is 99Dup plus the γ hydride `ZRH_B17` (a line compound). Hallstedt notes "other changes probably unintentional", for example the bcc Zr:H end member, which is based on GHSERZR instead of GBCCZR.
- **ZrH-07Uki** puts H in the tetrahedral sites: hcp (Zr)₁(H,Va)₂ and bcc (Zr)₁(H,Va)₆. The other Zr–H files use the octahedral sites, hcp 1:1 and bcc 1:3 (00Koe: bcc 1:1.5).
- **Laves phases** (NbZr-19Laf) are not stable in Nb–Zr and are not read: the reader takes only phases with mixing on one sublattice.

## Adding a database

Put the file here under the name it has on phasediagrams.org. Add a row to the table above. Check that the reader accepts the phases you need. It refuses anything outside its subset of the format with a message saying why.
