# Phototoxicity Datasets & In Silico Models — Sourced Survey

Compiled 2026-10-15. Focus: phototoxicity / photoirritation / photoallergy, with emphasis on the 3T3 Neutral Red Uptake (NRU) phototoxicity test (OECD TG 432) and related in vitro assays. Throughout, "open-source" = source code publicly available; "free to use" = usable at no cost without source access. Each fact carries its source URL, read date, and a `verified live` vs `index` flag.

## Summary

- **Phototoxicity data are scarce compared with Ames/MN.** The largest curated reference set found is the **PhotoChem Reference Chemical Database (251 compounds**, Toxics 2025), with in vitro (TG 432 / TG 495 / TG 498), in vivo, and human data. The regulatory anchor sets are tiny: **8 proficiency chemicals** in OECD TG 432 and **30 chemicals** in the ECVAM blind-trial validation study.
- **No standalone public photoallergy or photomutagenicity dataset** was found — both are gaps; photoallergy survives only as small literature-derived classification tables.
- **UV-Vis absorption data are the practical screening proxy** (ICH S10: MEC ≥ 1000 L/mol·cm in 290–700 nm), but absorption alone is necessary, not sufficient — a UV-Vis-based classifier predicted 3T3 NRU outcome with sensitivity 0.84 but specificity only 0.38.
- **Models are thin.** The only confirmed phototoxicity-specific predictors are **Derek Nexus (commercial, expert alerts)** and the free web tool **MLtox (QSAR + SHAP/XSMILES)**. **VEGA ships no phototoxicity model** (as far as its documented model set shows); **no phototoxicity-specific profiler was verified in the OECD QSAR Toolbox**; and **no open-source GitHub ML model trained specifically on phototoxicity data was found**.

---

## PART A — DATASETS

| Dataset | Endpoint(s) | ~Compounds | Format / access | Download URL | License / terms | Key assay details | Primary citation |
|---|---|---|---|---|---|---|---|
| OECD TG 432 proficiency chemicals | 3T3 NRU phototoxicity | 8 | Guideline table (not machine-readable) | OECD TG 432 PDF (see citation) | OECD | BALB/c 3T3 fibroblasts, NRU, ±UVA (5 J/cm²); PIF = IC50(−UV)/IC50(+UV); MPE = mean photo effect; thresholds: PIF <2 / MPE <0.1 non-phototoxic; PIF 2–5 / MPE 0.1–0.15 probable; PIF >5 / MPE >0.15 phototoxic | OECD TG 432 (2019) [index: https://www.oecd-ilibrary.org/test-no-432-in-vitro-3t3-nru-phototoxicity-test_5lmnfnj6bkq5.pdf] |
| ECVAM EU/COLIPA validation sets | 3T3 NRU phototoxicity | 30 (blind trial, 1992–1997) + 8 UV filters + 10 phototoxic + 10 non-phototoxic comparators (UV-filter special study) | Protocol PDF in ECVAM DB-ALM; results in cited papers (not a bulk download) | https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EURL-ECVAM/datasets/DBALM/LATEST/online/DBALM_docs/78_P_3T3%20Neutral%20Red%20Uptake%20(NRU)%20Phototoxicity%20Assay.pdf | Public (ECVAM/JRC) | Same 3T3 NRU PT protocol; blinded, multi-lab | Spielmann et al., Toxicol In Vitro 1998, 12:305–327 [index] |
| PhotoChem Reference Chemical Database | Phototoxicity + photoallergy (in vitro TG 432/495/498, in vivo, human) | 251 (human-verified: 59 phototoxic, 26 non-phototoxic, 4 photoallergenic; animal: 60 PT / 70 NPT) | Article tables (Table 3); standalone machine-readable SI NOT verified | https://www.mdpi.com/2305-6304/13/7/545 | Open access (MDPI CC BY) | CAS, physicochemical properties, UV absorption (λmax, ε), TG 432/495/498 results, human patch-test outcomes | Hong, Bae & Lim, Toxics 2025, 13(7):545 [verified live 2026-10-15] |
| MLtox training data | Phototoxicity (cell / 3D tissue / animal outcomes) | Size NOT verified | Not publicly downloadable (verified) | http://mltox.online (tool) | Not stated | Measured data on cell cultures, 3D models, animals | "MLtox, online phototoxicity prediction webpage", Toxicol In Vitro 2023 [index: https://www.sciencedirect.com/science/article/abs/pii/S0887233323001509] |
| UV-Vis / photoreactivity (Sci Rep 2021) | UV-Vis absorption features → photoreactive-potential class | 74,784 (train 72,788; test I 998; test II 998, incl. 43 with 3T3 NRU results) | Reaxys-derived; NO public download (proprietary source) | — | Reaxys terms (proprietary) | Methanol solvent; MW 98–1080; positive class = ≥1 absorption max 290–700 nm with MEC ≥ 1000 L/mol·cm (ICH S10 rule) | "Machine learning prediction of UV–Vis spectra features…", Sci Rep 2021 [index: https://www.nature.com/articles/s41598-021-03070-9] |
| UV/Vis spectra, experimental + computed (Sci Data 2019) | λmax + molar extinction coefficient (ε) | 8,488 parsed; 7,726 with experimental λmax (+ε); 762 computed-only | MongoDB/NoSQL hosting; bulk-download path NOT verified | — | Not verified | ChemDataExtractor-mined from literature; SMILES + solvents + units | "Comparative dataset of experimental and computational attributes of UV/vis absorption spectra", Sci Data 2019 [index: https://www.nature.com/articles/s41597-019-0306-0] |
| QM9S | Computed IR/Raman/UV-Vis spectra | 133,885 | figshare download | https://doi.org/10.6084/m9.figshare.24235333 | Public (figshare) | Computed spectra (DFT-based), NOT experimental | Nat Comput Sci 2023 (DetaNet/QM9S) [index] |
| ORNL_AISD-Ex / GDB-9-Ex | Computed excited states / UV-vis | Not verified (large) | Via Sci Data article | https://www.nature.com/articles/s41597-023-02408-4 | Public | Computed, NOT experimental | Sci Data 2023 [index] |
| ECHA REACH (registered substances / ECHA CHEM) | Phototoxicity endpoint summaries (OECD 432) | Per-substance (registration database) | https://echa.europa.eu/ (Information on Chemicals / ECHA CHEM) | Public factsheets; underlying study reports mostly confidential | Per-substance study summaries; no ECHA-published single bulk curated phototoxicity file found | REACH [index] |
| ChEMBL / PubChem bioassays | Phototoxicity | No curated phototoxicity assay ID verified | PubChem BioAssay / ChEMBL, queryable by assay type | https://pubchem.ncbi.nlm.nih.gov/ ; https://www.ebi.ac.uk/chembl/ | Public | Flag as gap: no specific phototoxicity AID verified — do not invent assay IDs | — |
| Photoallergy | Photoallergy | No standalone dataset found | Literature-derived classification tables only | See de Ávila et al. 2023, Cosmetics 11(2):47 (https://www.mdpi.com/2079-9284/11/2/47) | Open access | ~30+ compounds classified phototoxic+photoallergen / photoallergen; classification disagreements across studies noted | de Ávila et al. 2023 [index]; PhotoChem DB (4 human-verified photoallergenic) |
| Photomutagenicity | Photomutagenicity | No dataset found | — | — | — | Flag as gap | Lynch et al., Mutat Res 2011, 723:91–100 (context only) [index] |

Notes on specific datasets:
- **PhotoChem DB is the best starting point** for phototoxicity modeling: the largest curated set with paired in vitro + in vivo + human outcomes. But it is published as article tables, so reuse requires table extraction/curation (and a check of what the SI does/doesn't contain).
- **TG 432 proficiency chemicals and the ECVAM validation sets are reference/proficiency sets, not training sets** — far too small to train a model, but ideal as external validation / sanity checks.
- **The UV-Vis route is a proxy, not the endpoint.** The ICH S10 MEC ≥ 1000 rule and UV-Vis classifiers screen for photoreactivity; they do not measure phototoxicity. The Sci Rep 2021 test (n=43) shows the cost: sensitivity 0.84, specificity 0.38.
- **Reaxys-derived UV-Vis data cannot be redistributed** — do not treat the 74,784-compound set as an open dataset.

## PART B — PREDICTION MODELS / TOOLS

| Tool | Endpoints | Method | Training data | License | URL / repo | How to run |
|---|---|---|---|---|---|---|
| Derek Nexus (Lhasa) | **Phototoxicity** (+ phospholipidosis) among miscellaneous endpoints | Expert structural alerts + reasoning | Proprietary (Lhasa; incl. donated/pharma data) | **COMMERCIAL — not open** | https://www.lhasalimited.org/ [index] | Licensed desktop |
| MLtox | Phototoxicity | ML QSAR (molecular descriptors) + XSMILES/SHAP explainability | Measured cell / 3D-tissue / animal data (size not verified) | **FREE web server; NOT verified open-source** (no public code repo found) | http://mltox.online ; tool page https://mltox.fiit.stuba.sk/sk/data/predict/ [verified live 2026-10-15] | Web (input by CAS / SMILES / name); simple + advanced modes |
| VEGA (VEGAHUB) | **No phototoxicity model identified** in the documented model set | — | — | FREE to use; **not open-source** | https://www.vegahub.eu/ [verified live 2026-10-15] | GUI desktop app + online |
| OECD QSAR Toolbox | **No phototoxicity-specific profiler verified** | Read-across / profilers / databases (used alongside Derek Nexus in Ahuja et al. 2024) | Toolbox DBs | **FREE to download/use; NOT open-source** | https://qsartoolbox.org/ (v4.9) | Windows GUI; free registration/download |
| ProTox-3.0 | **No phototoxicity model identified** in the documented 61-model list | ML + similarity/fragment propensities | Curated in vitro/in vivo sets | **FREE web server (academic/non-commercial); NOT open-source** | https://tox.charite.de/protox3/ [index] | Web form only; free, no login |
| J Cheminform 2011 QSPR classifiers | Phototoxicity | PLS + recursive partitioning on quantum-chemical + classical descriptors (LUMO, HOMO/LUMO gap, electron affinity, logD, pKa…) | Phototoxicant/non-phototoxicant set (size not stated in abstract); external test set proprietary | Paper-only; **no public tool/code verified** | https://link.springer.com/article/10.1186/1758-2946-3-S1-O10 [index] | Not runnable as a tool |
| UV-Vis RF model (Sci Rep 2021) | UV-Vis photoreactive-potential class (proxy for phototoxicity screening) | Random forest on Morgan fingerprints | 72,788 train (Reaxys-derived) | Paper-only; **no public code/data verified** | https://www.nature.com/articles/s41598-021-03070-9 [index] | Not runnable as a tool |
| Open-source GitHub ML models trained on phototoxicity | — | — | — | — | — | **None found — gap.** General tox GNN work (e.g., ComptoxAI/qsar-gnn) covers Tox21 assays, not phototoxicity specifically |

## Caveats

1. **Tiny training data.** The largest curated phototoxicity set is 251 compounds; validation sets are 8–30. Any ML model trained here has a narrow applicability domain and high overfitting risk. Do not expect Ames-level performance.
2. **3T3 NRU overpredicts vs. in vivo.** In PhotoChem, TG 432 specificity against animal data was only 54.3% (19/35) — in vitro positives frequently do not reproduce in vivo. Models trained on 3T3 NRU labels inherit this bias.
3. **UV dose dependence.** Phototoxicity depends on the light dose (TG 432 uses 5 J/cm² UVA); models essentially never encode irradiance, exposure duration, or skin concentration. A "phototoxic" label is conditional on the assay's light dose.
4. **Solvent / assay variability.** Solvent choice measurably changes 3T3 NRU outcomes (e.g., bergamot oil: phototoxic in PBS, not in DMSO/ethanol). Labels pooled across solvents/protocols carry noise.
5. **Photoallergy has no validated in vitro method.** In silico photoallergy prediction is essentially unvalidated; do not present photoallergy classifications as model-validated.
6. **UV-Vis absorption is necessary, not sufficient.** MEC ≥ 1000 L/mol·cm in 290–700 nm is a screening rule (ICH S10). Absorption-positive compounds are often non-phototoxic (specificity 0.38 in the Sci Rep test) — absorption ≠ phototoxicity.
7. **Label disagreements across studies.** The photoallergy literature shows conflicting classifications for the same chemicals (e.g., acridine, 5-MOP, musk ketone) — curate calls per source and document resolution rules.
8. **Applicability domain.** Training sets skew toward pharmaceuticals/cosmetics ingredients (PhotoChem) or industrial chemicals (UV-Vis sets). Pharma-like, agrochemical, and nanomaterial space is unevenly covered; check AD before trusting a prediction.
9. **Derek Nexus is the only confirmed expert-alert system for phototoxicity, and it is commercial.** There is no open-source alert set for phototoxicity analogous to Toxtree's Benigni/Bossa rules for mutagenicity.

## Could not verify

- **PhotoChem DB machine-readable SI**: article tables confirmed; a downloadable SI file was not verified on the article page.
- **MLtox training-set size, model algorithm details, and code availability**: abstract + tool page only.
- **VEGA phototoxicity model**: none found in the documented model set — re-check in the VEGA app (model list evolves); absence is a finding from documentation, not a proof about the app binary.
- **QSAR Toolbox phototoxicity profiler**: none verified in documentation; the Ahuja et al. 2024 use of Toolbox was alongside Derek Nexus (databases/read-across), not a dedicated profiler.
- **ProTox-3.0 phototoxicity model**: none identified in the documented 61-model list.
- **ChEMBL / PubChem phototoxicity assay IDs**: none verified — do not use invented AIDs.
- **ECHA bulk phototoxicity extraction**: per-substance summaries are public; no ECHA-published single bulk curated phototoxicity file found.
- **Sci Data 2019 UV/Vis dataset bulk download**: MongoDB hosting described; direct download path not verified.
- **Live status of mltox.online (http)**: the STUBA tool page (mltox.fiit.stuba.sk) was verified live; the http://mltox.online domain from the paper was not separately opened.
- **Sarah Nexus / Leadscope / CASE Ultra phototoxicity coverage**: not verified — do not assume they ship phototoxicity models.

## Sources

All note files in `notes/` (one per source group, with URL, read date 2026-10-15, verified-live vs index flag): `datasets-3t3-photochem.md`, `datasets-uvvis-photoallergy-gap.md`, `models-phototox.md`. Verified live in-browser: https://www.mdpi.com/2305-6304/13/7/545 (PhotoChem article), https://www.vegahub.eu/ (home page), https://www.vegahub.eu/portfolio-item/vega-qsar/ (VEGA QSAR page), https://mltox.fiit.stuba.sk/sk/data/predict/ (MLtox tool page).
