# Data availability

This directory holds the small canonical inputs GloriPaka needs. The full provenance of
every dataset - what it measures, its publication and DOI, how it was harvested, cleaned,
deduplicated, scale-tagged, and rank-normalized - is documented in
[`docs/datasets.md`](../docs/datasets.md).

## What the package needs

GloriPaka predicts VHH HIC retention from sequence. It draws on two "worlds" of HIC data
that sit on non-comparable absolute scales (only within-assay rank is portable; see
`docs/datasets.md`, Section 1):

- **World 2 - the challenge scale (Merck ARRAS: TSKgel Butyl-NPR + ~50% methanol).**
  The 81 challenge VHH (`official_test_set.csv`, sequences only) and the native-VHH
  sources on or near this scale (Arras 2023; Siegmund 2025).
- **World 1 - salt-gradient HIC (Sepax Butyl-NP5 / PROPHET-Ab).** A large mAb/bsAb pool
  (Jain 2017, Shehata 2019, GDPa1/GDPa3/GDPa4) used for learning-curve and transfer
  analysis, plus the VHH sources Gordon 2026 / TNP and Ginkgo GDPa2.

Files this directory is expected to hold locally. None of them are committed to the
repository (everything under `data/` is git-ignored); obtain them from the sources below
and in `docs/datasets.md`:

| File | Contents |
|:--|:--|
| `official_test_set.csv` | The 81 challenge VHH sequences (no labels). |
| `native_vhh_labels.csv` | The 56-VHH native training pool: `Name, sequence, hic, hic_unit, scale_group, source_key, format`. |
| `hic_harvest/hic_records_vhh.csv` | Native-VHH provenance harvest with full assay context and `scale_group` tags. |
| `features/*.csv` | Feature tables computed by the project's own pipeline. |

## Third-party data is cited, not redistributed

Third-party HIC datasets remain the property of their originators and are governed by
their own licenses. This repository **cites** these datasets and records how to obtain
them; it does **not** redistribute third-party HIC measurements. Data the project
generated itself (predicted structures, derived feature tables) is not stored here either,
because it regenerates deterministically from sequence.

**Ginkgo / PROPHET-Ab data (GDPa1, GDPa2, GDPa3, GDPa4) is NOT included** in this
repository. To reproduce the World-1 mAb pool and the
GDPa2 portion of the native-VHH pool, obtain it directly from Ginkgo Datapoints. Any
GDPa2-derived rows present in a working copy of `hic_harvest/hic_records_vhh.csv` are
cite-only provenance records to be regenerated from the Ginkgo source, not a redistributed
release artifact.

## How to obtain each dataset

| Dataset | World | Type | Where to obtain | DOI |
|:--|:--:|:--|:--|:--|
| Arras 2023 (34 public labels) | 2 | VHH | *Front. Mol. Biosci.* article Table 4 (+ *mAbs* library paper) | 10.3389/fmolb.2023.1249247; 10.1080/19420862.2023.2261149 |
| Siegmund 2025 | 2 | VHH | *mAbs* SI (Fig S1 + Table S1) | 10.1080/19420862.2025.2553622 |
| Gordon 2026 / TNP | 1 | VHH | `github.com/oxpig/TNP` (repo CSVs) | 10.1038/s42003-026-09594-y |
| GDPa2 | 1 | VHH | `datapoints.ginkgo.bio` (**cite-only, not shipped**) | (verify) - Ginkgo portal / PROPHET-Ab platform |
| Jain 2017 | 1 | mAb | PNAS Datasets S2 + S3 | 10.1073/pnas.1616408114 |
| Shehata 2019 | 1 | mAb | Cell Reports SI (mmc2/3/5) | 10.1016/j.celrep.2019.08.056 |
| GDPa1 | 1 | mAb | `github.com/ginkgobioworks/abdev-benchmark` (**cite-only, not shipped**) | 10.1080/19420862.2025.2593055 |
| GDPa3 | 1 | mAb | `datapoints.ginkgo.bio` (**cite-only, not shipped**) | 10.1080/19420862.2026.2634216 (dataset: no separate DOI) |
| GDPa4 | 1 | bsAb | `datapoints.ginkgo.bio` (**cite-only, not shipped**) | (verify) |

DOIs marked **(verify)** could not be confirmed from the project's source records: GDPa2
has no confirmed dedicated publication or dataset DOI (it is a Ginkgo Datapoints portal
release under the PROPHET-Ab platform), and no GDPa4 publication DOI was confirmed.

For the complete dataset-by-dataset provenance, assay conditions, assembly methodology
(harvesting, sequence normalization, deduplication, the known/unknown split of the 81
test VHH, clone-family clustering, and within-`scale_group` rank normalization), and the
full summary table, see [`docs/datasets.md`](../docs/datasets.md).
