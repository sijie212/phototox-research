#!/usr/bin/env python3
"""Prepare the labeled phototoxicity dataset from the PhotoChem DB.

Source: Hong, Bae & Lim, "PhotoChem Reference Chemical Database for the
Development of New Alternative Photosafety Test Methods", Toxics 2025,
13(7):545, doi:10.3390/toxics13070545 — Table 3 (251 reference chemicals).

The article publishes the database as an article table (no machine-readable
SI; Data Availability: "available from the corresponding author on
request"). This script parses Table 3 from the open-access JATS XML via
Europe PMC (PMC12298655) and applies the label rule below.

LABEL RULE (documented, fixed):
  The endpoint modeled is the in vitro 3T3 Neutral Red Uptake (NRU)
  phototoxicity test (OECD TG 432), the standard in vitro phototoxicity
  assay. The TG432 column of Table 3 gives the call:
    TG432 == "PT"  -> label 1 (phototoxic)
    TG432 == "NPT" -> label 0 (non-phototoxic)
    anything else (blank, "PT/NPT" mixed, "-") -> EXCLUDED.
  Rationale: TG432 is the most complete column (149 clean calls vs. 85
  clean human calls and 130 clean animal calls), and it matches the
  assay-outcome-prediction framing used by the IVMN model. Human and
  animal columns are kept as provenance but do NOT define the label.

  Photoallergy: a distinct, immune-mediated endpoint (type IV
  hypersensitivity), NOT phototoxicity. This model predicts
  phototoxicity (TG432) only. Compounds annotated as photoallergenic in
  the Human/Animal columns enter under their TG432 call if they have a
  clean one (or are excluded if not); no separate photoallergy class is
  created. See MODEL.md for the count of such compounds.

SMILES: Table 3 has no SMILES column, so structures are resolved from
CAS numbers via PubChem PUG-REST (name lookup), cached in
data/cas_smiles_cache.json. CAS values that do not resolve, resolve to
an unparseable structure, or belong to mixtures/extracts/ill-defined
materials are dropped and counted. Output SMILES are RDKit-canonicalized
and deduplicated (first occurrence kept; conflicting-label duplicates
are reported and dropped).

Output: data/phototox_photochem.csv with columns:
  smiles, label, name, cas, human, animal, tg432, tg495, tg498
"""
import csv
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from xml.etree import ElementTree as ET

JATS_URL = ("https://www.ebi.ac.uk/europepmc/webservices/rest/"
            "PMC12298655/fullTextXML")
JATS_CACHE = "data/photochem_jats.xml"
CACHE = "data/cas_smiles_cache.json"
OUT = "data/phototox_photochem.csv"
TABLE_ID = "toxics-13-00545-t003"


def fetch_jats():
    if os.path.exists(JATS_CACHE):
        print(f"using cached {JATS_CACHE}")
        return open(JATS_CACHE, encoding="utf-8").read()
    print(f"downloading {JATS_URL}")
    req = urllib.request.Request(JATS_URL, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        text = resp.read().decode("utf-8")
    with open(JATS_CACHE, "w", encoding="utf-8") as fh:
        fh.write(text)
    return text


def parse_table3(xml_text):
    root = ET.fromstring(xml_text)
    for tw in root.iter("table-wrap"):
        if tw.get("id") == TABLE_ID:
            table = tw.find("table")
            rows = []
            for tr in table.findall(".//tr"):
                cells = ["".join(c.itertext()).strip() for c in tr]
                rows.append(cells)
            return rows
    raise RuntimeError(f"table {TABLE_ID} not found in JATS XML")


def resolve_cas_to_smiles(cas, cache):
    """PubChem PUG-REST name lookup by CAS. Returns SMILES or None."""
    if cas in cache:
        return cache[cas]
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + urllib.parse.quote(cas)
           + "/property/CanonicalSMILES/JSON")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        props = data["PropertyTable"]["Properties"]
        # NOTE (2026-10-15): PUG-REST currently returns the field as
        # "ConnectivitySMILES" (canonical connectivity SMILES, no stereo)
        # rather than "CanonicalSMILES". Accept either.
        smiles = (props[0].get("CanonicalSMILES")
                  or props[0].get("ConnectivitySMILES")
                  or props[0].get("IsomericSMILES")) if props else None
    except (urllib.error.HTTPError, urllib.error.URLError,
            KeyError, IndexError, ValueError, TimeoutError) as exc:
        print(f"  PubChem lookup failed for CAS {cas}: {exc}")
        smiles = None
    cache[cas] = smiles
    time.sleep(0.25)  # stay under PubChem's ~5 req/s guidance
    return smiles


def main():
    os.makedirs("data", exist_ok=True)
    rows = parse_table3(fetch_jats())
    header, body = rows[0], rows[1:]
    idx = {name: i for i, name in enumerate(header)}
    print(f"Table 3: {len(body)} data rows, columns: {header}")

    # ---- Label rule ----
    labeled = []
    n_missing = n_mixed = 0
    for r in body:
        tg432 = r[idx["TG432 (3T3)"]]
        if tg432 == "PT":
            label = 1
        elif tg432 == "NPT":
            label = 0
        else:
            if tg432.strip() == "":
                n_missing += 1
            else:
                n_mixed += 1
            continue
        labeled.append({
            "name": r[idx["Chemical"]],
            "cas": r[idx["CAS No."]],
            "human": r[idx["Human"]],
            "animal": r[idx["Animal"]],
            "tg432": tg432,
            "tg495": r[idx["TG495 (ROS)"]],
            "tg498": r[idx["TG498 (RhE Model)"]],
            "label": label,
        })
    n_pos = sum(r["label"] for r in labeled)
    print(f"TG432 label rule: {len(labeled)} labeled "
          f"({n_pos} PT / {len(labeled) - n_pos} NPT); "
          f"excluded {n_missing} missing + {n_mixed} mixed/other TG432")

    # Photoallergy annotation check (provenance only).
    pa = [r for r in labeled
          if "PS" in r["human"] or "PS" in r["animal"]]
    print(f"labeled rows with a photosensitivity (PS) annotation in "
          f"Human/Animal columns: {len(pa)} (kept under their TG432 call; "
          f"model predicts phototoxicity, not photoallergy)")

    # ---- CAS -> SMILES ----
    cache = {}
    if os.path.exists(CACHE):
        cache = json.load(open(CACHE))
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")

    seen = {}
    n_no_cas = n_unresolved = n_unparseable = n_dup = n_conflict = 0
    for r in labeled:
        cas = (r["cas"] or "").strip()
        if not cas or cas == "-":
            n_no_cas += 1
            continue
        raw = resolve_cas_to_smiles(cas, cache)
        if not raw:
            n_unresolved += 1
            continue
        mol = Chem.MolFromSmiles(raw)
        if mol is None:
            n_unparseable += 1
            continue
        smiles = Chem.MolToSmiles(mol)  # canonical
        if smiles in seen:
            n_dup += 1
            if seen[smiles]["label"] != r["label"]:
                n_conflict += 1
                print(f"WARN conflicting duplicate SMILES dropped: "
                      f"{smiles} (CAS {cas})")
                del seen[smiles]
            continue
        seen[smiles] = {"smiles": smiles, **r}
    with open(CACHE, "w") as fh:
        json.dump(cache, fh, indent=1, sort_keys=True)

    out_rows = sorted(seen.values(), key=lambda r: r["smiles"])
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["smiles", "label", "name", "cas",
                                           "human", "animal", "tg432",
                                           "tg495", "tg498"])
        w.writeheader()
        w.writerows(out_rows)
    n_pos = sum(r["label"] for r in out_rows)
    print(f"wrote {OUT}: {len(out_rows)} compounds "
          f"({n_pos} phototoxic / {len(out_rows) - n_pos} non-phototoxic)")
    print(f"dropped: {n_no_cas} no CAS, {n_unresolved} CAS unresolved, "
          f"{n_unparseable} unparseable SMILES, {n_dup} duplicate SMILES "
          f"({n_conflict} conflicting-label)")


if __name__ == "__main__":
    main()
