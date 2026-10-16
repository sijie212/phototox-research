# Phototoxicity Prediction — HF Embeddings + Few-Shot Learning

A binary classifier that predicts phototoxicity (`1` = phototoxic,
`0` = non-phototoxic) from a chemical structure (SMILES). Frozen
Hugging Face molecular embeddings (ChemBERTa) + a lightweight
classifier — no transformer fine-tuning. Built by directly reusing
the IVMN pipeline (same scripts, same CLI, same protocol, one
documented adaptation: the 128-sample few-shot size is infeasible at
this n — see Results).

## What it is

- **Endpoint:** in vitro 3T3 Neutral Red Uptake (NRU) phototoxicity
  test (OECD TG 432), binary call: `PT` (phototoxic) vs `NPT`
  (non-phototoxic).
- **Embeddings:** `DeepChem/ChemBERTa-77M-MLM` (RoBERTa, 384-dim),
  mean-pooled last hidden state, frozen. The embedding model is NOT in
  this repo: it downloads from Hugging Face on first run and is cached
  locally (`~/.cache/huggingface`) for offline reuse afterwards.
- **Classifier:** logistic regression (`C=1.0`, `class_weight='balanced'`)
  on the frozen embeddings. kNN (k=5, cosine) is evaluated alongside
  for comparison. The saved classifier is 2.5 KB (coefficients only).

## Dataset

**PhotoChem Reference Chemical Database**, Hong, Bae & Lim,
"PhotoChem Reference Chemical Database for the Development of New
Alternative Photosafety Test Methods", *Toxics* 2025, 13(7):545
(doi:10.3390/toxics13070545), Table 3 (251 reference chemicals).

The article publishes the database as an article table (no
machine-readable SI; Data Availability: "available from the
corresponding author on request"). `data/prepare_photochem.py` parses
Table 3 from the open-access JATS XML via Europe PMC (PMC12298655).

### Label rule (fixed, documented)

The TG432 column of Table 3 defines the label:

| TG432 call | Label | Handling |
|---|---|---|
| `PT` | 1 (phototoxic) | kept |
| `NPT` | 0 (non-phototoxic) | kept |
| blank / `PT/NPT` mixed / `-` | — | **excluded** |

Rationale: TG432 is the most complete column (149 clean calls vs. 85
clean human calls and 130 clean animal calls), and it matches the
assay-outcome-prediction framing of the IVMN model. The Human and
Animal columns are kept in the CSV as provenance but do **not**
define the label.

### Attrition (251 → 134)

| Step | Compounds |
|---|---|
| Table 3 rows | 251 (+ 6 header/section rows = 257 parsed rows) |
| Clean TG432 call (PT/NPT) | 149 (79 PT / 70 NPT) |
| Excluded: TG432 missing / mixed | 101 missing + 7 mixed/other |
| Excluded: no usable CAS | 1 |
| Excluded: CAS not resolvable in PubChem | 13 (mostly mixtures/extracts, e.g. bergamot oil, peppermint oil, "botanical drug") |
| Excluded: duplicate SMILES | 1 (same label on both rows) |
| **Final labeled set** | **134 (74 phototoxic / 60 non-phototoxic)** |

Table 3 has no SMILES column, so structures were resolved from CAS
numbers via PubChem PUG-REST (cached in
`data/cas_smiles_cache.json`), canonicalized with RDKit, and
deduplicated. One PUG-REST quirk worth knowing: the API currently
returns the field as `ConnectivitySMILES` rather than
`CanonicalSMILES`; `prepare_photochem.py` accepts either.

### Photoallergy is NOT what this model predicts

Photoallergy is a distinct, immune-mediated endpoint (type IV
hypersensitivity), not phototoxicity. This model predicts
**phototoxicity (TG432)** only. 12 of the 134 labeled rows carry a
photosensitivity (PS) annotation in the Human/Animal columns; they
enter under their TG432 call (or are excluded if they lack one) and
no separate photoallergy class is created.

## Install and run (fresh clone)

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt

# 1. Build the labeled dataset (parses PhotoChem Table 3 via Europe PMC,
#    resolves CAS -> SMILES via PubChem; cached afterwards)
.venv/bin/python data/prepare_photochem.py

# 2. Train + few-shot evaluation (embeds via ChemBERTa, ~1 min on CPU
#    after the one-time HF model download)
.venv/bin/python train_fewshot.py

# 3. Predict
.venv/bin/python predict.py "CCO" "Nc1ccc2ccccc2c1"
.venv/bin/python predict.py --csv input.csv --smiles-col smiles --out predictions.csv
```

First run needs internet (Hugging Face model download, ~300 MB, plus
the Europe PMC / PubChem data downloads). Behind a corporate proxy,
`embed.py` strips bracketed IPv6 entries from `NO_PROXY`/`no_proxy`
before importing `transformers`, because `huggingface_hub`'s HTTP
client crashes on them (`Invalid port: ':1]'`).

## Results

Fixed stratified split (seed 42): train 93 / val 20 / test 21.
Few-shot: 5 stratified subsets per size, drawn from the train split,
evaluated on the fixed test set. Mean ± std over the 5 repeats.

**Adaptation vs. the IVMN protocol:** the IVMN few-shot evaluation
also used size 128, but that is infeasible here — the phototox train
split (93 compounds) is smaller than 128, so a 128-sample subset
cannot be drawn. The final model, trained on train+val (113
compounds), is the max-data point instead.

### Few-shot test metrics (test n = 21)

| Train size | Classifier | Accuracy | ROC-AUC | Precision | Recall | F1 |
|---|---|---|---|---|---|---|
| 16 | logreg | 0.61 ± 0.07 | 0.70 ± 0.13 | 0.67 ± 0.07 | 0.67 ± 0.15 | 0.66 ± 0.07 |
| 16 | kNN | 0.61 ± 0.10 | 0.68 ± 0.12 | 0.63 ± 0.08 | 0.73 ± 0.22 | 0.67 ± 0.14 |
| 32 | logreg | 0.65 ± 0.13 | 0.75 ± 0.11 | 0.70 ± 0.11 | 0.63 ± 0.19 | 0.66 ± 0.15 |
| 32 | kNN | 0.66 ± 0.07 | 0.65 ± 0.07 | 0.66 ± 0.06 | 0.83 ± 0.11 | 0.73 ± 0.06 |
| 64 | logreg | 0.69 ± 0.05 | 0.71 ± 0.05 | 0.72 ± 0.04 | 0.73 ± 0.06 | 0.73 ± 0.05 |
| 64 | kNN | 0.60 ± 0.04 | 0.62 ± 0.06 | 0.63 ± 0.05 | 0.77 ± 0.06 | 0.69 ± 0.02 |

### Final model (logreg trained on train+val, 113 compounds; test n = 21)

| Model | Accuracy | ROC-AUC | Precision | Recall | F1 |
|---|---|---|---|---|---|
| **logreg (shipped)** | **0.71** | **0.77** | **0.75** | **0.75** | **0.75** |
| kNN (comparison) | 0.52 | 0.64 | 0.57 | 0.67 | 0.62 |

Confusion matrix (logreg): [[TN=6, FP=3], [FN=3, TP=9]].

**Read these numbers with the sample size in mind.** Test n = 21
means one compound flipping changes accuracy by ~5 points, and the
few-shot standard deviations (±0.05–0.15 AUC) are wide. The headline
AUC of 0.77 is encouraging for such a small dataset, but it is a
noisy estimate — not a precise measurement. See Caveats.

## Interpretability (MCTS)

`interpret_mcts.py` finds a **rationale**: a minimal substructure of
an input molecule that preserves the model's prediction — the
fragment the prediction "rests on". It is a deletion-based Monte
Carlo Tree Search, same implementation as the IVMN model:

- **State:** a SMILES string, starting from the full input molecule.
- **Action:** delete one heavy atom and keep the largest connected
  fragment (RDKit). Empty/invalid fragments are skipped; fragments
  whose aromaticity breaks on deletion are repaired via Kekulization.
- **Reward:** the model's predicted probability of the target class
  for the fragment, minus a small size penalty (`--size-penalty`,
  default 0.05), so smaller fragments that keep the prediction win.
- **Search:** MCTS with UCB1 selection (`--simulations`, default
  200). Nodes with at most `--max-atoms` heavy atoms (default 12)
  are terminal. Each new fragment costs one ChemBERTa forward pass,
  so a molecule finishes in ~1–2 min on CPU.
- **Reported rationale:** the smallest evaluated fragment that
  preserves the predicted class, verified to be a valid SMILES and a
  substructure of the input. If no fragment preserves the prediction,
  the script says so instead of returning a bogus rationale.

```bash
.venv/bin/python interpret_mcts.py "CC(=O)Oc1ccccc1C(=O)O"
.venv/bin/python interpret_mcts.py "CCO" --target 0 --simulations 100 \
    --out rationale.json
```

Output: input prediction + probability, the rationale SMILES with
its heavy-atom count and probability, and the top candidate
fragments by reward (with whether each keeps the prediction).
`--out` writes the same result as JSON. `--target {0,1}` overrides
which class to explain (default: the model's predicted class).

Verified examples (200 simulations each, independently re-checked
for valid SMILES + substructure match + class preservation):

| Molecule | Model call | Rationale | Rationale P(pos) |
|---|---|---|---|
| Protoporphyrin IX (label 1) | positive, 0.70 | 30-atom porphyrin core fragment (input 42) | 0.65 |
| Aspirin (label 0) | negative, 0.35 | `O=CO` (3 atoms) | 0.45 |
| Ethanol | negative, 0.21 | `CCO` (whole molecule, 3 atoms) | 0.21 |

The input probabilities above are the model's calls on the full
molecules as reported by the script runs.

**The Aspirin rationale is degenerate — and that is informative.**
A 3-atom formic-acid fragment "preserving" the negative call does
not mean formic acid drives the prediction; it means the model's
negative call on aspirin is not anchored in any specific moiety —
nearly any small fragment stays on the negative side. Rationale
quality tracks model quality: with a test AUC of 0.77 on n = 21,
rationales are sanity-check devices, not mechanistic evidence.

**Limitation — read this before citing a rationale.** A rationale is
a *model-explanation* device: it shows which part of the structure
the model relies on to keep its prediction. It is **not** proof of a
biological mechanism — it does not show that the fragment causes
phototoxicity.

## Caveats — read before using

- **Tiny data, wide uncertainty.** 134 training compounds; test
  n = 21. Every metric carries a wide confidence interval. One
  borderline compound flipping changes accuracy by ~5 points. Do
  not quote AUC 0.77 as a precise figure.
- **Label-definition sensitivity.** The label is the TG432 (3T3 NRU)
  call. TG432 specificity against animal data was only 54.3% in the
  PhotoChem analysis — in vitro positives frequently do not
  reproduce in vivo. A model trained on TG432 labels inherits this
  bias. Human and animal labels would give a different (smaller)
  dataset and possibly different rankings; this choice is documented
  so it can be revisited.
- **Light dose is not encoded.** Phototoxicity depends on the light
  dose (TG 432 uses 5 J/cm² UVA); the model sees only structure,
  never irradiance, exposure duration, or skin concentration. A
  "phototoxic" label is conditional on the assay's light dose.
- **Solvent / assay variability.** Solvent choice measurably changes
  3T3 NRU outcomes (e.g., bergamot oil: phototoxic in PBS, not in
  DMSO/ethanol). Labels pooled across solvents/protocols carry
  noise. 13 CAS entries (mostly mixtures/extracts) could not be
  resolved to a single structure and were dropped.
- **Probabilities near 0.5 are low-confidence.** Treat borderline
  calls as coin flips, not real predictions.
- **Applicability domain.** 134 mostly pharmaceutical/cosmetic
  compounds; industrial, agrochemical, and nanomaterial space is
  unevenly covered. Check similarity to the training set before
  trusting a prediction.
- **Screening baseline, not a regulatory replacement.** For ICH S10
  photosafety decisions, follow the established testing battery.
  Photoallergy is not covered by this model at all.

## Files

| File | What it is |
|---|---|
| `embed.py` | HF embedding extraction (ChemBERTa, mean-pooled) |
| `train_fewshot.py` | Few-shot evaluation + final model training; writes `metrics.json` |
| `predict.py` | CLI / CSV prediction with the saved classifier |
| `interpret_mcts.py` | MCTS rationale extraction (minimal class-preserving substructure) |
| `model_phototox.joblib` | Saved logreg classifier (2.5 KB) |
| `metrics.json` | Few-shot + final test metrics (machine-readable) |
| `MODEL.md` | This file (named MODEL.md because the repo root already has a survey README.md) |
| `requirements.txt` | Pinned deps (torch CPU, transformers, scikit-learn, rdkit, …) |
| `data/prepare_photochem.py` | Parses PhotoChem Table 3 → `data/phototox_photochem.csv` |
| `data/photochem_jats.xml` | Cached PhotoChem JATS XML from Europe PMC (source of Table 3) |
| `data/cas_smiles_cache.json` | Cached PubChem CAS → SMILES resolutions |
| `data/phototox_photochem.csv` | Labeled dataset (134 compounds) |
| `data/splits.csv` | Train/val/test assignment (seed 42) |
| `data/embeddings.npz` | Cached ChemBERTa embeddings (134 × 384, regenerable via `embed.py`) |
