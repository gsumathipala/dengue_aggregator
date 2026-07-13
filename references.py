"""
Dengue Genome Aggregator - Reference Genomes and Protein Coordinates

Loads the 4 DENV serotype RefSeq reference genomes (GenBank flat files under
references/) and parses their mature-peptide (mat_peptide) feature tables to
get real, per-serotype protein boundaries.

NOTE: the reference accessions here were independently verified against NCBI
Taxonomy (2026-07-12). The DENV_REFERENCES dict previously hardcoded in
variant_annotator.py was wrong for 3 of 4 serotypes:
  - serotype 1 was assigned NC_001474, which is actually DENV-2's RefSeq.
  - serotype 2 was assigned NC_001438, which is a Bean golden yellow mosaic
    virus DNA-B segment (unrelated plant virus) -- not dengue at all.
  - serotype 3 was assigned NC_001469, which is Cucumber necrosis virus --
    also not dengue.
  - serotype 4's NC_002640 was correct.
Confirmed correct accessions (via NCBI taxon IDs 11053/11060/11069/11070,
cross-checked against esummary titles):
    DENV-1 -> NC_001477.1 (10735 bp)
    DENV-2 -> NC_001474.2 (10723 bp)
    DENV-3 -> NC_001475.2 (10707 bp)
    DENV-4 -> NC_002640.1 (10649 bp)
"""

import logging
import re
from pathlib import Path
from typing import Dict, Optional, Tuple

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

REFERENCES_DIR = Path(__file__).parent / "references"

REFERENCE_ACCESSIONS = {
    1: "NC_001477",
    2: "NC_001474",
    3: "NC_001475",
    4: "NC_002640",
}

# Canonical mature-peptide product names (as they appear verbatim in the
# RefSeq GenBank records) mapped to short protein labels. Both DENV
# polyprotein precursor fragments ("anchored capsid protein ancC", the
# standalone "protein pr", and standalone "membrane glycoprotein M") are
# intentionally excluded in favor of the mature "capsid protein C" and the
# aggregate "membrane glycoprotein precursor prM", matching the 10(+1)-region
# scheme used elsewhere in this codebase (C, prM, E, NS1, NS2A, NS2B, NS3,
# NS4A, 2K, NS4B, NS5).
PRODUCT_NAME_TO_PROTEIN = {
    "capsid protein c": "C",
    "membrane glycoprotein precursor prm": "prM",
    "envelope protein e": "E",
    "nonstructural protein ns1": "NS1",
    "nonstructural protein ns2a": "NS2A",
    "nonstructural protein ns2b": "NS2B",
    "nonstructural protein ns3": "NS3",
    "nonstructural protein ns4a": "NS4A",
    "protein 2k": "2K",
    "nonstructural protein ns4b": "NS4B",
    "rna-dependent rna polymerase ns5": "NS5",
}

_cache: Dict[int, dict] = {}


def _parse_genbank_sequence(text: str) -> str:
    """Extract the nucleotide sequence from a GenBank flat file's ORIGIN block."""
    start = text.index("ORIGIN")
    end = text.index("\n//", start)
    body = text[start:end]

    bases = []
    for line in body.split("\n"):
        m = re.match(r"^\s*\d+\s+([ACGTNacgtn \t]+)$", line)
        if m:
            bases.append(re.sub(r"[^ACGTNacgtn]", "", m.group(1)))
    return "".join(bases).upper()


def _parse_genbank_mat_peptides(text: str) -> Dict[str, Tuple[int, int]]:
    """
    Extract mat_peptide (start, end) coordinates keyed by short protein name.

    Returns 0-indexed, Python-slice-style (start inclusive, end exclusive)
    coordinates converted from GenBank's 1-indexed inclusive coordinates.
    """
    features_start = text.index("FEATURES")
    features_end = text.index("ORIGIN")
    section = text[features_start:features_end]

    proteins = {}
    for m in re.finditer(
        r"mat_peptide\s+(\d+)\.\.(\d+)[\s\S]*?/product=\"([^\"]+)\"", section
    ):
        genbank_start, genbank_end, product = m.group(1), m.group(2), m.group(3)
        protein = PRODUCT_NAME_TO_PROTEIN.get(product.strip().lower())
        if protein is None:
            continue
        start_0idx = int(genbank_start) - 1
        end_0idx_exclusive = int(genbank_end)
        proteins[protein] = (start_0idx, end_0idx_exclusive)

    return proteins


def get_reference(serotype: int) -> dict:
    """
    Load (and cache) the reference genome and protein map for a DENV serotype.

    Returns:
        {
            'accession': str,
            'sequence': str,
            'length': int,
            'proteins': {protein_name: {'start': int, 'end': int, 'length': int}},
        }
    """
    if serotype in _cache:
        return _cache[serotype]

    accession = REFERENCE_ACCESSIONS.get(serotype)
    if not accession:
        raise ValueError(f"No reference defined for serotype {serotype}")

    path = REFERENCES_DIR / f"{accession}.gb"
    if not path.exists():
        raise FileNotFoundError(
            f"Reference GenBank file not found: {path}. "
            f"Fetch it with efetch (db=nuccore, id={accession}, rettype=gb)."
        )

    text = path.read_text(encoding="utf-8")
    sequence = _parse_genbank_sequence(text)
    mat_peptides = _parse_genbank_mat_peptides(text)

    proteins = {
        name: {"start": start, "end": end, "length": end - start}
        for name, (start, end) in mat_peptides.items()
    }

    missing = set(PRODUCT_NAME_TO_PROTEIN.values()) - set(proteins)
    if missing:
        logger.warning(f"DENV-{serotype} ({accession}): missing expected protein regions: {missing}")

    result = {
        "accession": accession,
        "sequence": sequence,
        "length": len(sequence),
        "proteins": proteins,
    }
    _cache[serotype] = result
    logger.info(f"Loaded DENV-{serotype} reference {accession} ({len(sequence)} bp, "
                f"{len(proteins)} annotated protein regions)")
    return result


def get_protein_at_position(serotype: int, position: int) -> Tuple[Optional[str], Optional[int]]:
    """
    Determine which mature protein a 0-indexed reference genome position falls into.

    Returns:
        (protein_name, protein_start) or (None, None) if intergenic/UTR.
    """
    ref = get_reference(serotype)
    for protein_name, coords in ref["proteins"].items():
        if coords["start"] <= position < coords["end"]:
            return protein_name, coords["start"]
    return None, None


if __name__ == "__main__":
    for st in (1, 2, 3, 4):
        ref = get_reference(st)
        print(f"DENV-{st}: {ref['accession']}, {ref['length']} bp, proteins={sorted(ref['proteins'])}")
