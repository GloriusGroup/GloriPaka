# Datasets and data provenance

This document records, for every hydrophobic-interaction-chromatography (HIC) dataset
used in the GloriPaka project, where the datapoints came from, what they measure, the
source publication and DOI, how the records were obtained and assembled (harvesting,
sequence cleaning, deduplication, scale tagging, rank-normalization), and their
redistribution status.

Everything below is grounded in the project's own provenance records and per-source
literature notes. Where a DOI, count, or figure could not be confirmed from those
sources it is marked **(verify)** rather than guessed. No third-party datapoints
(sequences paired with HIC values) are copied into this document; each dataset is
described, cited, and pointed to its origin.

**Source records used to compile this document**

- `modeling_strategyB/notes/lab_journal.md` - the primary running record: counts, scale
  groups, assembly decisions, and the known/unknown split.
- `modeling_strategyB/README.md` - the "two-worlds scale reality", dataset counts, and
  scale bridges.
- `data/SOURCES.md` - download dates, URLs, and licenses (all downloads dated 2026-07-01).
- Assembled tables: `data/hic_harvest/hic_records_vhh.csv`,
  `modeling_strategyB/data/native_vhh_labels.csv`,
  `modeling_strategyB/data/world1_labels.csv` / `world1_input.csv`,
  `data/HIC_datasets/ginkgo_gdpa/scale_bridges.csv`, `data/official_test_set.csv`.
- Per-source literature notes in `literature/notes/` (Arras 2023, Jain 2017,
  Arsiwala 2025 / PROPHET-Ab, van Niekerk 2025/2026, Gordon 2026 / TNP, Siegmund 2025,
  Waibl 2022, and the data-mining / pretrain-pool inventories).

A one-page data-availability statement lives in `data/README.md`.

---

## 1. Overview: the two-scale reality

HIC retention time is an assay readout, not an intrinsic molecular constant. It depends
on the column chemistry, the salt system, the presence or absence of an organic
modifier, the gradient length, flow rate, temperature, detection wavelength, and the
analyte format. As a direct consequence, **absolute HIC retention times are not
comparable across different assays; only the rank order within a single run/assay is
portable.** The clinical-antibody landscape paper (Jain 2017) already reports only
Spearman rank correlations between assays for exactly this reason, and the 2025 Ginkgo
Datapoints competition was designed around a single platform specifically to remove
inter-laboratory assay heterogeneity as a confounder (van Niekerk 2025/2026).

The project therefore separates all HIC data into two "worlds":

- **World 1 - salt-gradient HIC (PROPHET-Ab / Sepax-Butyl style, no organic modifier).**
  A large pool of conventional monoclonal antibodies (mAb) and bispecific antibodies
  (bsAb): Jain 2017, Shehata 2019, and the Ginkgo PROPHET-Ab series (GDPa1, GDPa3,
  GDPa4). These sit on the Sepax Proteomix Butyl-NP5 column family (ammonium-sulfate
  gradient to phosphate buffer). They are large and internally rank-comparable but are
  **not** the challenge scale. World 1 is used for learning-curve and transfer analysis,
  not for the headline VHH model.

- **World 2 - Merck ARRAS scale (TSKgel Butyl-NPR + ~50% methanol).** The challenge
  target scale: a Tosoh TSKgel Butyl-NPR column eluted with an ammonium-sulfate to
  50%-methanol gradient. Only two sources live natively on this scale: Arras 2023
  (the source of the 34 public test labels) and Siegmund 2025 (three VHH-Fc constructs
  on a byte-for-byte identical assay). There is **no numeric bridge** from World 1 to
  World 2; transfer is by rank and by biophysical features only.

### Why absolute minutes are not comparable, and how the project handles it

Because absolute minutes are not portable, the project never pools raw retention times
across assays. Instead it:

1. Tags every record with a `scale_group` string identifying the exact assay it was run
   on (see the table in Section 2).
2. Rank / quantile-normalizes **within** each `scale_group` before any pooling
   (`within_group_quantiles` in the modeling code). This maps each source onto a common
   rank axis so that VHH measured on three different machines can train one model
   without any source's absolute scale leaking into another (lab journal, Entries 10 and
   16).
3. Uses **scale bridges** - sets of molecules measured on two scales at once - to
   confirm cross-scale rank equivalence where it exists. The bridges are stored in
   `data/HIC_datasets/ginkgo_gdpa/scale_bridges.csv` (columns:
   `bridge, molecule, scale_a, value_a, scale_b, value_b`; 202 paired rows):

   | Bridge | Molecules (n) | Scales linked | Spearman rho |
   |:--|:--:|:--|:--:|
   | `Jain2017 <-> GDPa1_B` | 133 | Sepax Butyl-NP5 20-min (Jain, scale A) vs 4.5-min (GDPa1, scale B) | **0.968** |
   | `GDPa4_norm <-> GDPa1_B` | 69 | PROPHET-Ab normalized RT vs PROPHET-Ab raw RT | **0.986** |

   These high Spearman values show the Sepax-Butyl-NP5 gradients are essentially one
   rank axis (a time-compressed version of one another), which is why World-1 mAb
   sources pool cleanly as a rank block. **No such bridge exists to the ARRAS/methanol
   World-2 scale** - the two worlds share no molecules measured on both, so World 1 to
   World 2 transfer is done through sequence-derived biophysical features, not through
   any numeric calibration.

---

## 2. Assay-scale groups (as tagged in the assembled tables)

The `scale_group` column is the canonical assay tag. The values that actually occur in
the assembled files are:

| scale_group | World | Column / eluent | Gradient | Sources tagged with it |
|:--|:--:|:--|:--|:--|
| `ARRAS_TSKButylNPR_MeOH` | 2 | Tosoh TSKgel Butyl-NPR; 1.2 M (NH4)2SO4/PBS to 50% methanol/0.1x PBS; UV 214 nm | 15 min, 0.5 mL/min, 25 C | Arras 2023 (37), Siegmund 2025 (3) |
| `SEPAX_ButylNP5_salt` | 1 | Sepax Proteomix (HTC) Butyl-NP5; 1.8 M (NH4)2SO4/0.1 M phosphate pH 6.5 to phosphate; UV 280 nm | 20 min, 1.0 mL/min | Gordon 2026 / TNP (36) |
| `GINKGO_GDPa2` | 1 | Ginkgo PROPHET-Ab HIC-HPLC (Sepax Butyl-NP5 family); compressed 0-4.5-min scale | 4.5-min class | GDPa2 (51 rows / 17 VHH) |
| `SEPAX_ButylNP5_gradientA` | 1 | Sepax Proteomix Butyl-NP5; ~20-min salt gradient pH 6.5, UV 280 nm | ~20 min, 1 mL/min | Jain 2017 (137), Shehata 2019 (400) |
| `SEPAX_ButylNP5_gradientB` | 1 | Sepax Proteomix Butyl-NP5; 4.5-min PROPHET-Ab gradient | 4.5 min, 1.5 mL/min | GDPa3 (80); GDPa1 (`GINKGO_PROPHET_B`) is the same gradient |
| `GINKGO_PROPHET_B` / `_norm` | 1 | PROPHET-Ab raw / normalized RT | 4.5-min class | GDPa1 (242), GDPa4 (160 bsAb + 71 mAb re-measures) |

Note: Siegmund's three VHH-Fc constructs are tagged to the ARRAS scale group because
their HIC assay is byte-for-byte identical to the Arras assay (same column, pre-mix,
gradient, standards). GDPa2 is kept on its own `GINKGO_GDPa2` group because it is a
single-replicate PROPHET-Ab run and mixes two analyte formats (bare-VHH-His and
VHH-Fc), so it is not assumed equal to the mAb PROPHET-Ab scales until proven.

---

## 3. Datasets, one subsection each

### 3.1 Official challenge test set - 81 VHH (World-2 target)

- **What it is.** The 81 VHH sequences the project must predict, `data/official_test_set.csv`
  (columns: `ID`, `VHHsequence`; 81 rows, ~80 unique sequences). Sequences only; no HIC
  labels are provided by the challenge.
- **What it measures.** Nothing directly - these are the query molecules. Their true HIC
  values are measured on the World-2 ARRAS assay (see 3.2).
- **Molecule type.** Single-domain antibodies (VHH / nanobody), Merck's humanized
  anti-NKp46 platform, assayed as VHH-Fc fusions. Framework analysis (lab journal,
  Entry 3) confirms all 81 number cleanly IMGT 1 to 128 (they are full V-domains, not
  fragments); 71 of the 81 carry a VHH FR2 hallmark tetrad, 10 carry a conventional
  VH-type FR2 (`DVQL...`, humanized/universal-framework).
- **Provenance.** Strong evidence (lab journal Entry 1/3; `literature/notes/arras_2023.md`)
  that these are Merck Healthcare KGaA (Evers/Zielonka) VHH from the same platform as
  Arras 2023: the shared framework (`EVQLLESGGG...WGQGTLVTVSS` and the `DVQLQESGGG...`
  sibling germline) matches the Arras constructs exactly.
- **Known / unknown split.** Exact sequence matching of the 81 against the native-VHH
  harvest yields **34 "known"** (all on `ARRAS_TSKButylNPR_MeOH`, HIC 4.7-6.7 min) and
  **47 "unknown"**. The 34 are largely parent clones; the 47 include engineered
  point-mutant variants of them. The 34 are used strictly as a held-out validation proxy
  and never enter training, feature selection, or normalization statistics (lab journal
  "Cardinal rule"). A looser Strategy-A CDR-combination overlap counts up to 44 rows
  (`data/official_test_hic_*` files); those are leakage files and are not used here.
- **Publication / DOI.** Challenge-provided; provenance publications are the two Arras
  2023 papers (3.2). The test file itself has no DOI.
- **Redistribution.** Challenge input; shipped as sequences only (no labels).

### 3.2 ARRAS - Arras et al. 2023 (World-2 native scale; source of the 34 public labels)

Two companion papers from Merck Healthcare KGaA / EMD Serono (Evers & Zielonka group)
describing one project: a llama immunized against human NKp46, CDR3 repertoire grafted
onto a humanized VHH scaffold, yeast-display selection, then one-armed VHH-SEEDbody
Fc-fusion production and developability profiling including analytical HIC.

- **Paper 1 (source of `arras37`).** Arras P, Yoo HB, Pekar L, Clarke T, Friedrich L,
  Schröter C, Schanz J, Tonillo J, Siegmund V, Doerner A, Krah S, Guarnera E, Zielonka S,
  Evers A. "AI/ML combined with next-generation sequencing of VHH immune repertoires
  enables the rapid identification of de novo humanized and sequence-optimized single
  domain antibodies: a prospective case study." *Front. Mol. Biosci.* 10:1249247 (28 Sep
  2023). **DOI: 10.3389/fmolb.2023.1249247.**
- **Paper 2 (library approach; further ~42 VHH not yet mined into the pool).** Arras P,
  Yoo HB, Pekar L, Schröter C, Clarke T, Krah S, Klewinghaus D, Siegmund V, Evers A,
  Zielonka S. "A library approach for the de novo high-throughput isolation of humanized
  VHH domains with favorable developability properties following camelid immunization."
  *mAbs* 15(1):2261149 (27 Sep 2023). **DOI: 10.1080/19420862.2023.2261149** (Correction:
  10.1080/19420862.2023.2285575). Open Access CC BY-NC.
- **What it measures / N.** Paper 1 Table 4 gives **37 VHH with HIC** (IDs 1-37, range
  4.7-6.7 min); these are the `arras2023_frontiers` rows in the harvest (37 records).
  Paper 2 Table 3 gives a further **42 VHH with HIC** (range 4.5-7.3 min) that are largely
  non-redundant but **not yet incorporated** into the assembled pool (they require SI
  Table S4 for full WT/VGLW sequences). Together the two papers hold ~79 VHH on the ARRAS
  assay.
- **Molecule type.** VHH, assayed as one-armed VHH-SEEDbody Fc fusions.
- **Assay (World-2 target scale).** Tosoh TSKgel Butyl-NPR, 2.5 µm, 4.6 x 100 mm; 20 µg
  pre-mixed 1:1 with 2 M ammonium sulfate; linear gradient A = 1.2 M ammonium sulfate in
  PBS to B = 50% methanol in 0.1x PBS, 15 min, 25 C, 0.5 mL/min; UV 214 nm; standards
  cetuximab (5.7-5.8 min) and avelumab (7.2 min).
- **How obtained.** The 37 HIC values were transcribed from Paper 1 Table 4; full amino
  acid sequences were reconstructed from the fixed humanized FERF framework plus the
  Table 3 CDR1/CDR2/CDR3 definitions (values cross-checked against `arras37.csv`, e.g.
  ID1 = 4.9, ID15 = 6.4, ID31 = 6.4, ID36 = 6.7). Sequence-to-HIC pairing was verified in
  `literature/notes/arras_reconciliation_verified.md`.
- **Role.** The 34 of these that exactly match the test set are the held-out validation
  proxy; the full 37 are featurized separately (`arras37_all.csv`). Arras is the
  best-matched in-scale anchor for the challenge.
- **Redistribution.** Cite-only. Values live in the article tables/SI; obtain from the
  DOIs above. No public raw-HIC repository exists beyond the article.

### 3.3 Siegmund et al. 2025 (World-2 native scale; 3 VHH-Fc in the training pool)

- **Publication.** Siegmund V, Klewinghaus D, Teroerde J, Pekar L, Harwardt J, Zielonka S,
  Militano F, Licari G, Evers A. "Optimizing colloidal stability and viscosity of
  multispecific antibodies at the drug discovery-development interface: a systematic
  predictive case study." *mAbs* 17(1):2553622 (2025). **DOI: 10.1080/19420862.2025.2553622.**
  Open Access CC BY-NC. Merck Healthcare KGaA (same group and same HIC assay as Arras).
- **What it measures / N.** HIC (and colloidal-stability) measurements on 15 engineered
  constructs (9 bispecific IgG-VHH + 3 Fc-VHH + 3 Fab-Fc IgG controls). Only **3 VHH-Fc**
  constructs enter the native-VHH training pool (`siegmund2025`, format
  "VHH-Fc (N-term Fc fusion)"), HIC 4.4-7.2 min.
- **Molecule type.** VHH domains in Fc-VHH / bispecific IgG-VHH format.
- **Assay.** Byte-for-byte identical to Arras 2023 (TSKgel Butyl-NPR; 2 M ammonium
  sulfate 1:1 pre-mix; A = 1.2 M ammonium sulfate/PBS to B = 50% methanol/0.1x PBS,
  15 min; avelumab anchored at 7.2 min). This is one of the very few external HIC sets
  that is directly on the challenge scale.
- **How obtained.** Sequences from Figure S1 and HIC retention times from Table S1 of the
  SI, paired manually (the paper deposits no data in external repositories).
- **Caveat.** All 15 are pI-engineering mutants of one parent, so sequence diversity is
  near-zero; value is a few on-scale points plus assay confirmation, not a training-set
  expander. The paper's mechanistic finding - VHH HIC is driven by charge distribution
  and conformational stability, not hydrophobicity alone - informed the charge/pI
  features.
- **Redistribution.** Cite-only. Obtain sequences and HIC values from the mAbs SI at the
  DOI above.

### 3.4 Gordon et al. 2026 - Therapeutic Nanobody Profiler (TNP) (native VHH, World-1 scale)

- **Publication.** Gordon GL, Gervasio J, Souders C, Deane CM. "Characterising nanobody
  developability to improve therapeutic design using the Therapeutic Nanobody Profiler."
  *Communications Biology* 9:344 (2026). **DOI: 10.1038/s42003-026-09594-y.** Open Access
  CC BY 4.0. Oxford OPIG + Twist Bioscience. Code/data: `github.com/oxpig/TNP` (repo
  license BSD-3).
- **What it measures / N.** 108 nanobodies with experimental HIC retention time (105
  valid), of which **36 are clinical-stage VHH with public sequences** and 72 are
  proprietary Twist VHH with HIC but withheld sequences. Only the **36 clinical VHH**
  enter the training pool (`gordon2026_tnp`, HIC 11.8-20.9 min in the assembled table).
- **Molecule type.** VHH, assayed as VHH-IgG1 Fc fusions.
- **Assay (World-1 salt-gradient).** Sepax Proteomix HTC Butyl-NP5, 4.6 x 100 mm; 5 µg
  diluted to ~1 M ammonium sulfate; A = 1.8 M ammonium sulfate/0.1 M phosphate pH 6.5 to
  B = 0.1 M phosphate pH 6.5 (no organic modifier); linear 20 min, 1.0 mL/min, UV 280 nm.
  Retention 11.8-27 min - an entirely different absolute scale from ARRAS (~4.5-7.3 min);
  usable only after per-assay rank normalization.
- **How obtained.** Two repo CSVs were joined (verified all 36 match, 0 missing):
  in-vitro HIC from `paper/paper_data/invitro_descriptors/108_Nbs_in_vitro_data.csv`
  (`HIC Retention Time (min)` column) joined on drug name (`Samples` == `SeqID`) to
  sequences in `paper/paper_data/insilico_descriptors/VHH_TSD_all_properties_FINAL.csv`.
- **Redistribution.** Cite-only. Obtain from the oxpig/TNP GitHub repository; the 72
  proprietary sequences are unrecoverable (withheld).

### 3.5 GDPa2 - Ginkgo PROPHET-Ab VHH release (native VHH, own scale)

- **Publication / DOI.** Ginkgo Datapoints portal release (catalog title "18 VHH
  Constructs"), part of the PROPHET-Ab series. The platform is described by Arsiwala et al.
  2025 (*mAbs* 17(1):2593055, DOI 10.1080/19420862.2025.2593055), **but that paper
  contains no VHH data** - GDPa2 has no confirmed dedicated publication or dataset DOI.
  **DOI: (verify)** - treat GDPa2 as a Ginkgo Datapoints portal release under the
  PROPHET-Ab platform. (An earlier data-mining note recorded "GDPa2 not published"; it was
  subsequently released on the portal, so that note is superseded.)
- **What it measures / N.** 51 measured construct-variant rows corresponding to **17
  distinct parent nanobodies** (structure: 17 parents x 3 rows = 1 His-tagged bare-VHH +
  2 identical-sequence VHH-Fc replicate rows). HIC range 0.68-3.21 min on the compressed
  PROPHET-Ab scale. The training pool uses the **17 bare-VHH rows** (`ginkgo_gdpa2`,
  format "bare-VHH", HIC 0.68-2.84 min in the assembled table).
- **Molecule type.** VHH (clinical/benchmark nanobodies), measured as monovalent
  single-domain constructs (bare-VHH-His and VHH-Fc); several clinical parents are
  multivalent in their true format, but GDPa2 measured the single-domain building block.
- **Assay (own scale).** Ginkgo PROPHET-Ab HIC-HPLC (Sepax Butyl-NP5 family, compressed
  0-4.5-min axis). Single technical replicate (n=1) per value, so noisier than
  GDPa1/3/4. VHH-Fc RT > bare-VHH RT for the same VHH, so the two formats are **not**
  pooled on one axis.
- **How obtained.** Standardized into the project's HIC schema and appended to the
  native-VHH harvest (`hic_records_vhh.csv`, rows R0077-R0127, `scale_group=GINKGO_GDPa2`).
- **Redistribution.** **Cite-only. NOT redistributed here.** Obtain from Ginkgo Datapoints
  (`datapoints.ginkgo.bio` / Hugging Face `ginkgo-datapoints`). See Section 6.

### 3.6 Jain et al. 2017 (World-1 mAb; PNAS clinical-stage landscape)

- **Publication.** Jain T, Sun T, Durand S, Hall A, Houston NR, Nett JH, Sharkey B,
  Bobrowicz B, Caffry I, Yu Y, Cao Y, Lynaugh H, Brown M, Baruah H, Gray LT, Krauland EM,
  Xu Y, Vásquez M, Wittrup KD. "Biophysical properties of the clinical-stage antibody
  landscape." *PNAS* 114(5):944-949 (31 Jan 2017). **DOI: 10.1073/pnas.1616408114.**
  All authors employed by Adimab LLC (single-lab study).
- **What it measures / N.** 137 clinical-stage mAbs with in-house-measured HIC retention
  time (assembled as `jain`, 137 rows; 135 usable, 2 right-censored). Sequences were
  data-mined (WHO-INN 92 / US patent 26 / PDB 19); all biophysical values including HIC
  are fresh, standardized, single-lab Adimab measurements (not compiled from other
  papers).
- **Molecule type.** Full IgG mAb (paired VH/VL), re-expressed as human IgG1.
- **Assay (World-1, scale A).** Sepax Proteomix HIC Butyl-NP5; A = 1.8 M ammonium sulfate
  + 0.1 M sodium phosphate pH 6.5 to B = 0.1 M sodium phosphate pH 6.5; ~1 M ammonium
  sulfate at load; linear gradient over 20 min, 1 mL/min, UV 280 nm; 5 µg IgG; RT = major
  peak. Method taken from Estep et al. 2015 (see 3.11).
- **Censoring.** A retention of **25.00 min is a right-censored "no elution"** flag
  (Dataset S3 footnote a: "Arbitrarily long RT of 25 min to indicate no elution"), not a
  measured value. Such rows are treated as censored, not regressed on 25.00.
- **How obtained.** HIC values from PNAS Dataset S3 (`pnas.1616408114.sd03.xlsx`),
  sequences from Dataset S2 (`pnas.1616408114.sd02.xlsx`). (The Hugging Face mirror
  `hugging-science/jain-clinical-antibody-polyreactivity` retains only a binary
  polyreactivity label with HIC stripped, so it is not used for HIC.)
- **Redistribution.** Cite-only. Obtain the two PNAS supplementary xlsx files from the
  article at the DOI above.

### 3.7 Shehata et al. 2019 (World-1 mAb; repertoire panel)

- **Publication.** Shehata L, et al. (repertoire mAb biophysical study), *Cell Reports*
  (2019). **DOI: 10.1016/j.celrep.2019.08.056** (PMID 31553901). CC BY 4.0 (per Europe PMC).
  Full author list and article title **(verify)** - not transcribed in the project notes.
- **What it measures / N.** A large repertoire panel of mAbs with HIC retention time. The
  assembled World-1 label file (`world1_labels.csv`) carries **400** Shehata rows
  (`source=shehata`); at source level, Data S1 has ~402 repertoire mAbs of which ~348 have
  a valid HIC value (54 not determined), plus Data S2 (42 clinical anchors, ~36 usable)
  and Data S4 / a VH1-69 engineered panel (26 rows, ~18 usable) available in the SI but
  only partially merged. The lab-journal cross-panel analysis uses "Shehata 2019 (348)"
  as the valid-HIC repertoire count.
- **Molecule type.** Repertoire mAb (paired VH/VL); entirely novel sequences (Jain and
  Shehata share zero exact VH strings; the 42 clinical anchors share drug names but are
  re-expressed, so their sequence strings differ).
- **Assay (World-1, scale A).** Sepax Proteomix Butyl-NP5, ~20-min salt gradient pH 6.5,
  UV 280 nm - the same column family and gradient class as Jain (numerically similar,
  8.5-13.9 min). Same separation mechanism, different lab, so still rank-pooled within its
  own scale group.
- **How obtained.** Downloaded from the Cell Reports supplementary xlsx (mmc2 = Data S1,
  mmc3 = Data S2, mmc5 = Data S4); assembled locally to
  `data/HIC_datasets/shehata2019_mab_hic.csv`. (The Hugging Face mirror
  `hugging-science/shehata-antibody-psr` is polyreactivity-only, no HIC.)
- **Redistribution.** Cite-only. Obtain the supplementary xlsx from the Cell Reports
  article at the DOI above.

### 3.8 GDPa1 - Ginkgo PROPHET-Ab (World-1 mAb; primary mAb training set)

- **Publication.** Arsiwala A, Bhatt R, van Niekerk L, Quintero-Cadena P, Ao X,
  Rosenbaum A, ... Tessier PM, Borhani DW. "A high-throughput platform for biophysical
  antibody developability assessment to enable AI/ML model training." *mAbs* 17(1):2593055
  (02 Dec 2025). **DOI: 10.1080/19420862.2025.2593055.** Open Access CC BY-NC (paper).
  Platform = PROPHET-Ab. Code: `github.com/ginkgobioworks/prophet_ab`.
- **Dataset license.** The GDPa1 dataset is released **CC BY 4.0 with restrictions**
  (commercial use permitted except selling/licensing the dataset itself); see `LICENSE.md`
  in the Hugging Face repo.
- **What it measures / N.** 246 clinical/approved IgGs assayed on PROPHET-Ab; **HIC
  populated for 242** (this is the "242 mAb with HIC" figure). Range 2.43-4.50 min, mean
  2.82; only 3 mAbs at the 4.5-min cap (assigned RT = gradient end, analogous to Jain's
  censoring). Assembled as `gdpa1_all.csv` (242 mAb feature table); the raw CSV is
  `GDPa1_v1.2_20250814.csv` (v1.3 re-release also referenced).
- **Molecule type.** 100% full IgG (IgG1 169 / IgG2 31 / IgG4 46; kappa 217 / lambda 29).
  No VHH content.
- **Assay (World-1, scale B).** Sepax Proteomix Butyl-NP5 (identical column chemistry and
  salt system to Jain); A = 1.8 M ammonium sulfate/0.1 M phosphate pH 6.5 to B = 0.1 M
  phosphate pH 6.5; linear 4.5-min gradient, 1.5 mL/min, UV 280 nm; 2.5 µg full IgG.
  Test-retest reliability Pearson r = 0.99 (N = 237). GDPa1 HIC is effectively a
  time-compressed Jain HIC (rank rho = 0.97 vs Jain, n = 133).
- **How obtained.** Downloaded 2026-07-01 from the non-gated GitHub
  `ginkgobioworks/abdev-benchmark` (raw base
  `raw.githubusercontent.com/ginkgobioworks/abdev-benchmark/main/`); the HF mirror
  `ginkgo-datapoints/GDPa1` requires login and was not used. The CSV ships VH/VL/HC/LC
  sequences, AHo alignments, and three CV-fold columns (random, hierarchical-cluster, and
  the recommended isotype-stratified clustered fold).
- **Redistribution.** **Cite-only. NOT redistributed here.** Obtain from the Ginkgo
  source (GitHub `abdev-benchmark` or `datapoints.ginkgo.bio`). See Section 6.

### 3.9 GDPa3 - Ginkgo blinded competition hold-out (World-1 mAb)

- **Publication (context).** van Niekerk L, Moller J, Ritter S, et al. "2025 Ginkgo
  Datapoints Antibody Developability Competition outcomes: limited model performance and a
  call for data standardization." *mAbs* 18(1):2634216 (22 Feb 2026).
  **DOI: 10.1080/19420862.2026.2634216** (Correction: 10.1080/19420862.2026.2644040). Open
  Access CC BY-NC. GDPa3 is the competition hold-out; it has no separate dataset DOI.
- **What it measures / N.** 80 blinded IgGs (79 with HIC), OAS-sampled for diversity and
  disjoint from GDPa1, assayed on the same PROPHET-Ab platform; released publicly after the
  competition. HIC 2.34-4.80 min (3 at the 4.8 cap). Assembled as `gdpa3` (80 rows) in
  `world1_labels.csv`, `scale_group=SEPAX_ButylNP5_gradientB`.
- **Molecule type.** Repertoire mAb (paired VH/VL).
- **Assay.** PROPHET-Ab (scale B), same as GDPa1; note the cap differs (4.8 vs 4.5), so
  the two are the same column family but not byte-identical scales (cap values kept
  censored).
- **How obtained.** Downloaded from `datapoints.ginkgo.bio/dataset-access` (free).
- **Redistribution.** **Cite-only. NOT redistributed here.** Obtain from Ginkgo Datapoints.

### 3.10 GDPa4 - Ginkgo PROPHET-Ab bispecific + mAb re-measures (World-1 bsAb)

- **Publication / DOI.** Ginkgo Datapoints PROPHET-Ab release. A publication DOI of the
  form `10.1080/19420862.2026.*` is referenced in the project notes but not confirmed -
  **DOI: (verify)**.
- **What it measures / N.** 160 bispecific IgG (bsAb) with HIC (`GDPa4.1 bsAb`,
  `hihplc_normretentiontime_median`, 2.42-4.80 min) plus **71** monospecific IgG re-measures
  (69 with HIC) that are VH+VL-exact duplicates of GDPa1 molecules - i.e. 0 new molecules;
  they exist to bridge scales (see the GDPa4_norm <-> GDPa1_B bridge, Section 1).
- **Molecule type.** Bispecific IgG (bsAb, all four chains present) + duplicate mAb.
- **Assay.** PROPHET-Ab **normalized** HIC RT (`GINKGO_PROPHET_B_norm`), its own sub-scale
  (hence the `_norm` tag). The 160 bsAb are a new molecule class and are not rank-poolable
  with mAb until anchored via the 69-mAb overlap set.
- **Role.** Not in the headline VHH model; used as a scale bridge and a distinct bsAb class.
- **Redistribution.** **Cite-only. NOT redistributed here.** Obtain from Ginkgo Datapoints.

### 3.11 Supporting / candidate sources (not in the assembled pool)

- **Estep et al. 2015** - *mAbs* 7(3):553-561, DOI **(verify)**. The HIC method source
  cited by Jain 2017 (also an Adimab paper). ~32 approved therapeutics on a Dionex ProPac
  HIC-10 column (a different column, scale group C), sequence-redundant with Jain ->
  ~0 new molecules. Not in the pool. Cite-only.
- **Rauschmayer et al. 2026** - "A side-by-side biophysical comparison of five
  single-domain antibody scaffolds used for synthetic display libraries," *J. Biol. Chem.*
  302(7):113155, **DOI 10.1016/j.jbc.2026.113155** (PMC13284441). 5 sdAb with HIC (Thermo
  mAbPac HIC-10, its own scale); sequences require SI/origin papers. Flagged as a small
  candidate (`candidates/`) but not merged into the training pool. Cite-only.
- **Waibl et al. 2022** - "Comparison of hydrophobicity scales for predicting biophysical
  properties of antibodies," *Front. Mol. Biosci.* 9:960194, **DOI 10.3389/fmolb.2022.960194.**
  A methods/benchmark paper on hydrophobicity scales (used to select the HIC-tuned scales
  in the feature stack), not a training dataset. Its Roche-34/Roche-127 HIC data are
  partially proprietary and not used.
- **A dedicated 2026 re-scout** (`literature/notes/more_vhh_hic_data_candidates.md`)
  confirmed that no genuinely new, public, mineable VHH+HIC sequence-to-value dataset
  exists beyond Arras, Gordon/TNP, Siegmund, GDPa2, and Rauschmayer. The public VHH+HIC
  well is, to a good approximation, dry; the ceiling is set by data availability on a
  shared assay, not by method.

---

## 4. Assembly methodology (cross-cutting)

### 4.1 Harvesting workflow

All downloads were performed 2026-07-01 via `curl` from public, non-gated URLs where
possible; gated sources (e.g. the HF GDPa1 mirror) were avoided in favor of non-gated
equivalents (the `abdev-benchmark` GitHub). Every source is standardized into a common
HIC schema. The native-VHH harvest `data/hic_harvest/hic_records_vhh.csv` (127 rows) is
the canonical record; its columns are:

`record_id, source_key, citation, doi_or_patent, molecule_id, molecule_class, format,
sequence_region, sequence, hic_value, hic_unit, readout_type, censored, assay_column,
assay_salt, assay_modifier, assay_gradient_min, assay_pH, assay_detect_nm, ref_standards,
scale_group, provenance, notes`.

Each row therefore carries not just a sequence and HIC value but the full assay context
(`assay_column` ... `assay_detect_nm`), a censoring flag, a citation + DOI, the analyte
`format`, and the `scale_group` tag that drives all downstream normalization. The
harvest currently contains: Arras 2023 (37), Gordon/TNP (36), Siegmund 2025 (3), and
GDPa2 (51 rows / 17 VHH).

### 4.2 Sequence normalization

Sequences are stored as full V-domains. IMGT numbering (via ANARCI / anarcii) confirmed
that all VHH sources (test set, ARRAS, GDPa2, TNP, Siegmund) number cleanly IMGT 1 to
128, so no padding or trimming of the test set is needed (lab journal, Entry 3). Where a
source ships Fc-fusion constructs (for example GDPa2 VHH-Fc rows, ~357-369 aa), the Fc is
stripped by keeping only the ANARCI V-domain span before featurization; the raw fusion is
never fed to the structure or embedding models. For GDPa2 the bare-VHH-His rows were used
directly so no Fc-strip was required.

### 4.3 Deduplication policy

Deduplication is by **exact V-domain sequence string** plus normalized molecule name:

- GDPa1 is a superset of Jain: 135 of Jain's 137 VH strings are byte-identical in GDPa1;
  Jain contributes only 2 unique molecules over GDPa1, and GDPa1 adds 111 new clinical
  mAbs.
- Jain and Shehata share zero exact VH strings; the 42 Shehata clinical anchors share drug
  names but were re-expressed, so their sequence strings differ and they are treated as
  distinct.
- Estep's ~32 approved therapeutics are sequence-redundant with Jain (net ~0 new).
- GDPa4's 71 monospecific IgG are VH+VL-exact re-measures of GDPa1 (0 new molecules; kept
  only as a scale bridge).
- The 34 "known" test VHH are matched to the harvest by exact sequence; engineered
  point-mutant variants deliberately fall on the "unknown" side so they cannot leak a
  near-label.

### 4.4 How the native-VHH training pool was built

The native-VHH training pool (`modeling_strategyB/data/native_vhh_labels.csv`; columns
`Name, sequence, hic, hic_unit, scale_group, source_key, format`) is **56 unique VHH**
across three public sources and three assay scales:

| source_key | n | format | scale_group | HIC range (min) |
|:--|:--:|:--|:--|:--:|
| `gordon2026_tnp` | 36 | VHH-IgG1 Fc fusion | `SEPAX_ButylNP5_salt` | 11.8-20.9 |
| `ginkgo_gdpa2` | 17 | bare-VHH | `GINKGO_GDPa2` | 0.68-2.84 |
| `siegmund2025` | 3 | VHH-Fc (N-term Fc fusion) | `ARRAS_TSKButylNPR_MeOH` | 4.4-7.2 |

Arras 2023's 37 VHH are deliberately **excluded** from this training pool (they are the
validation proxy and are featurized separately as `arras37_all.csv`), so the pool is
disjoint by source from the 34 held-out labels - enabling train-once/predict-the-34 with
zero leakage (lab journal, Entry 10). Enlarging the pool from an early 34-VHH/6-family
version to this 56-VHH/3-source version is what flipped native-VHH training from a
negative to a positive result.

The World-1 mAb pool is assembled analogously in `world1_labels.csv` /
`world1_input.csv` (columns `Name, sequence, hic, source, scale_group` and
`Name, Heavy_Chain`): 617 deduped VH domains = GDPa3 (80) + Shehata (400) + Jain (137),
each tagged with its scale group (`SEPAX_ButylNP5_gradientA` for Jain/Shehata,
`SEPAX_ButylNP5_gradientB` for GDPa3). GDPa1 (242) is featurized separately, bringing the
World-1 mAb total to roughly 860 measurement rows (~756 unique molecules); GDPa4's 160
bsAb are a separate bispecific class.

### 4.5 Known/unknown split of the 81 test VHH

Defined by **exact V-domain sequence match** of the 81 against the native-VHH harvest:
matches on the ARRAS scale become the 34 "known" (HIC 4.7-6.7 min), the remaining 47 are
"unknown". This split is fixed in one place (`labels.known_test_labels()`) and reused
everywhere. The 34 known ARRAS labels are the project's held-out ruler: they never enter
training, feature selection, normalization statistics, or hyper-parameter choice.

### 4.6 Clone-family clustering for leak-free cross-validation

Because some of the hidden 47 are engineered near-duplicates of the known 34, honest CV
requires that near-duplicate clones never straddle a train/validation fold. Sequences are
grouped into clone families by pairwise identity. A naive 0.90 identity threshold merges
distinct clones (two different VHH already share ~100 framework residues, so cross-clone
identity is ~0.85-0.92); calibrating against the nearest-neighbor identity distribution of
the 81 (min 0.731, median 0.968) led to a chosen threshold of **0.95**, which merges true
variants (differing by a few residues) while keeping distinct clones apart (41 families
across the 81; 6 families across the 34 known). Leave-clone-family-out CV then uses these
groups (lab journal, Entry 1).

### 4.7 Rank-normalization used to pool across scales

Pooling across scales is always done by rank within a `scale_group`
(`within_group_quantiles`): each source's HIC values are mapped onto the training
target's empirical quantiles (a monotone transform, hence Spearman-invariant) before any
model sees them. This is the single most important data-handling rule in the project: it
lets VHH measured on three different machines train one model, and it guarantees that no
source's absolute minutes can leak a scale into another. Where a numeric cross-scale link
does exist within World 1, it is documented as a scale bridge (Section 1) rather than used
to rescale raw values into the model.

---

## 5. Data-availability and redistribution statement

**Third-party datasets are cited, not redistributed.** GloriPaka builds on public
datasets that remain the property of their originators and are governed by their own
licenses. This repository and any associated Zenodo deposit do not redistribute
third-party HIC measurements; they record provenance and give instructions to obtain each
dataset from its original source. Only data that the project generated itself and that is
license-clear (for example predicted structures and derived feature tables computed by the
project's own pipeline) may be shared via Zenodo.

**Ginkgo / PROPHET-Ab data specifically (GDPa1, GDPa2, GDPa3, GDPa4).** The Ginkgo
Datapoints / PROPHET-Ab measurements are **NOT included in this repository or in any
associated Zenodo deposit.** To reproduce the World-1 mAb pool and the GDPa2 portion of
the native-VHH pool, obtain them directly from Ginkgo:

- GDPa1: `github.com/ginkgobioworks/abdev-benchmark` (non-gated) or
  `datapoints.ginkgo.bio` / Hugging Face `ginkgo-datapoints/GDPa1` (login). Dataset license
  CC BY 4.0 with restrictions.
- GDPa2 (VHH), GDPa3, GDPa4: `datapoints.ginkgo.bio/dataset-access` and the
  `ginkgo-datapoints` Hugging Face organization.

Any development copy of `data/hic_harvest/hic_records_vhh.csv` that carries GDPa2-derived
rows is a working provenance record only; those rows are cite-only Ginkgo measurements and
are excluded from the redistributable release. Regenerate them from the Ginkgo source.

**Other third-party sources.** Jain 2017 (PNAS SI), Shehata 2019 (Cell Reports SI), Gordon
2026 / TNP (`github.com/oxpig/TNP`), Arras 2023 (Frontiers + mAbs article tables/SI), and
Siegmund 2025 (mAbs SI) are likewise cite-only: fetch them from the DOIs/repositories in
Section 3. The 81 challenge test sequences are challenge input (sequences only, no labels).

---

## 6. Summary table

Redistribution status: "cite-only" = obtain from the original source (not shipped here);
"input" = the sequences-only challenge file. Molecule types: mAb, bsAb, VHH. "N" is the
count as recorded in the project's assembled files / notes.

| Dataset | World | Molecule | N (HIC) | Assay conditions | Publication (DOI) | How obtained | Redistributable here |
|:--|:--:|:--|:--:|:--|:--|:--|:--|
| Official test set (81 VHH) | 2 | VHH | 81 seq (34 labels public) | target = ARRAS TSKgel Butyl-NPR + 50% MeOH | Merck challenge; provenance Arras 2023 | provided by challenge | input (seqs only) |
| Arras 2023 (Frontiers) | 2 | VHH | 37 | TSKgel Butyl-NPR; 1.2 M AS/PBS to 50% MeOH; 15 min; UV 214 | Arras 2023, *Front. Mol. Biosci.* 10:1249247 (10.3389/fmolb.2023.1249247) | Table 4 + framework reconstruction | cite-only |
| Arras 2023 (library, mAbs) | 2 | VHH | 42 (not yet merged) | same as above | Arras 2023, *mAbs* 15(1):2261149 (10.1080/19420862.2023.2261149) | Table 3 + SI Table S4 | cite-only |
| Siegmund 2025 | 2 | VHH | 3 (of 15) | identical to Arras assay | Siegmund 2025, *mAbs* 17(1):2553622 (10.1080/19420862.2025.2553622) | SI Fig S1 (seq) + Table S1 (HIC) | cite-only |
| Gordon 2026 / TNP | 1 | VHH | 36 (of 108) | Sepax HTC Butyl-NP5; 1.8 M AS to phosphate pH 6.5; 20 min; UV 280 | Gordon 2026, *Commun. Biol.* 9:344 (10.1038/s42003-026-09594-y) | join of two `oxpig/TNP` repo CSVs | cite-only |
| GDPa2 (Ginkgo) | 1 | VHH | 17 (51 rows) | PROPHET-Ab HIC (Sepax Butyl-NP5 family), 4.5-min class, n=1 | Ginkgo portal / PROPHET-Ab platform (Arsiwala 2025); own DOI (verify) | Ginkgo Datapoints portal | cite-only (NOT shipped) |
| Jain 2017 | 1 | mAb | 137 (135 usable) | Sepax Butyl-NP5; 1.8 M AS/phosphate pH 6.5; 20 min; UV 280; 25 min = censored | Jain 2017, *PNAS* 114(5):944 (10.1073/pnas.1616408114) | PNAS Datasets S2 (seq) + S3 (HIC) | cite-only |
| Shehata 2019 | 1 | mAb | 400 (~348 valid) | Sepax Butyl-NP5; ~20-min salt gradient pH 6.5; UV 280 | Shehata 2019, *Cell Reports* (10.1016/j.celrep.2019.08.056; PMID 31553901); title/authors (verify) | Cell Reports SI xlsx (mmc2/3/5) | cite-only |
| GDPa1 (Ginkgo) | 1 | mAb | 242 (of 246) | PROPHET-Ab: Sepax Butyl-NP5; 1.8 M AS/phosphate pH 6.5; 4.5-min; UV 280 | Arsiwala 2025, *mAbs* 17(1):2593055 (10.1080/19420862.2025.2593055) | GitHub `abdev-benchmark` (non-gated) | cite-only (NOT shipped) |
| GDPa3 (Ginkgo) | 1 | mAb | 80 (79 HIC) | PROPHET-Ab (scale B), 4.8-min cap | van Niekerk 2025/2026, *mAbs* 18(1):2634216 (10.1080/19420862.2026.2634216); dataset no separate DOI | `datapoints.ginkgo.bio` | cite-only (NOT shipped) |
| GDPa4 (Ginkgo) | 1 | bsAb (+mAb dup) | 160 bsAb (+71 mAb) | PROPHET-Ab normalized HIC RT (`_norm`) | Ginkgo PROPHET-Ab release; DOI (verify) | `datapoints.ginkgo.bio` | cite-only (NOT shipped) |
| Estep 2015 | (C) | mAb | ~32 (~0 new) | Dionex ProPac HIC-10 (different column) | Estep 2015, *mAbs* 7(3):553-561; DOI (verify) | article (method source for Jain) | cite-only |
| Rauschmayer 2026 | (other) | VHH/sdAb | 5 (candidate) | Thermo mAbPac HIC-10 | Rauschmayer 2026, *JBC* 302(7):113155 (10.1016/j.jbc.2026.113155) | article + SI (candidate, not merged) | cite-only |

Data the project generated itself (predicted NanoBodyBuilder2 structures, computed
feature tables, rank-normalization outputs) is license-clear and may be shared via Zenodo;
it is derived from, but does not reproduce, any third-party HIC measurement.
