# GloriPaka — packaging roadmap

Turning the Strategy-B research workspace (`modeling_strategyB/`) into a clean,
reproducible, standalone Python package. Checked items are done; the rest are the
plan. Built **minimal-first, then expand**.

## Phase 1 — Minimal viable repo (installable + importable) ✅ in progress
- [x] Clean repo tree under `gloripaka/`; package renamed `strategyb` → `gloripaka`.
- [x] `config.py` decoupled from the parent Merck repo (self-contained, repo-relative paths).
- [x] `pyproject.toml` for **uv** (metadata, deps, optional extras, ruff+mypy+pytest, Google docstring convention).
- [x] `LICENSE` (MIT), `CITATION.cff`, `.gitignore`, `.gitattributes`, `py.typed`.
- [x] `README.md` with logo, badges, TOC, install (uv), quick start, module table, citation.
- [ ] `uv lock` + `uv sync` green; `python -c "import gloripaka"` works.
- [ ] `git init`, remote → `GloriPaka.git`, first commit; **push** (auth via personal PAT/SSH — deploy keys disabled).

## Phase 2 — End-to-end pipeline (sequence → prediction)
- [ ] Port featurization (`featurization/`) into `src/gloripaka/featurize/`; dedup vs `featurization_w1/`; strip `__pycache__`; keep DeepSP vendored weights + license.
- [ ] `src/gloripaka/structure.py` — NanoBodyBuilder2 structure generation from sequence.
- [ ] `src/gloripaka/pipeline.py` — high-level `GloriPakaPredictor`: **sequences in → structures → features → HIC + confidence out**.
- [ ] CLI entry point (`gloripaka predict input.csv`).
- [ ] Curate `experiments/` (keep the E0–E23 keepers; each runs against in-repo data).
- [ ] `scripts/reproduce.py` — regenerates the headline result + key figures.

## Phase 3 — Data + Zenodo
- [ ] Consolidate the small canonical CSVs into `data/` (inputs, labels, feature tables actually used).
- [ ] Keep a tiny in-repo sample (34 known VHH) under `data/examples/` for tests/examples.
- [ ] Stage the heavy artifacts (PDBs, full feature tables) in the separate `zenodo_gloripaka/` folder; write deposition metadata; upload; add DOI + `download_data.py`.
- [ ] `data/README.md` documenting provenance + how to fetch.

## Phase 4 — Literature + dataset construction (the "why")
- [ ] Curate `literature/notes/` learnings into `docs/literature.md` (cited).
- [ ] `docs/datasets.md` — how each dataset (World-1 pool, native VHH, scale bridges, ARRAS test) was constructed, with DOIs.
- [ ] Decide manuscript inclusion **carefully** (journal/figure license) — likely link out, not commit.

## Phase 5 — Quality: tests, typing, docs, CI
- [ ] `tests/` — smoke + unit (metrics, physics composite, within-group quantiles, clone-family clusterer, model fit/predict on toy).
- [ ] Docstring/typing pass → `ruff check` + `mypy src/` clean.
- [ ] Docs site (mkdocs-material) → GitHub Pages (docs badge).
- [ ] GitHub Actions: ruff + mypy + pytest on push/PR (tests badge).
- [ ] `CONTRIBUTING.md`, issue templates.

## Phase 6 — Release
- [ ] Tag `v0.1.0`, GitHub release; finalize Zenodo/docs/CI badges; README figures; cross-links.

---
_Source of truth for the science: `modeling_strategyB/notes/lab_journal.md`. This
roadmap tracks only the packaging/reproducibility effort._
