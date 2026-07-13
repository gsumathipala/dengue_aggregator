"""
Dengue Genome Aggregator - Variant Annotation Engine
Aligns each genome to its serotype-specific reference genome and calls
SNPs/indels with protein, codon, and functional-impact context.

Alignment uses rapidfuzz's bit-parallel global edit-distance implementation
(Levenshtein.editops) rather than a naive same-length position comparison --
real GenBank/BV-BRC/Nextstrain sequences are not pre-aligned to the reference
and vary in length, so comparing them position-by-position without aligning
first would misattribute every base downstream of any indel.

Reference genomes and their real, per-serotype protein coordinates (parsed
from RefSeq GenBank annotations) live in references.py.
"""

import logging
from dataclasses import dataclass
from typing import List, Optional

from rapidfuzz.distance import Levenshtein

import references

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


@dataclass
class Variant:
    """Variant record (reference-genome coordinates, 1-indexed)."""
    position: int
    ref_allele: str
    alt_allele: str
    variant_type: str  # SNP or INDEL
    impact: str  # HIGH, MODERATE, LOW, MODIFIER
    protein: Optional[str]
    protein_position: Optional[int]
    codon_change: Optional[str]
    aa_change: Optional[str]
    is_synonymous: bool
    is_frameshift: bool
    is_stop_gained: bool
    is_stop_lost: bool

    def to_dict(self) -> dict:
        return {
            'position': self.position,
            'ref_allele': self.ref_allele,
            'alt_allele': self.alt_allele,
            'variant_type': self.variant_type,
            'impact': self.impact,
            'protein': self.protein,
            'protein_position': self.protein_position,
            'codon_change': self.codon_change,
            'aa_change': self.aa_change,
            'is_synonymous': self.is_synonymous,
            'is_frameshift': self.is_frameshift,
            'is_stop_gained': self.is_stop_gained,
            'is_stop_lost': self.is_stop_lost
        }


class VariantAnnotator:
    """Aligns a query genome to its serotype's reference and calls variants."""

    def __init__(self, serotype: int):
        """
        Args:
            serotype: DENV serotype (1-4). The reference genome and protein
                      coordinates are loaded automatically from references.py.
        """
        self.serotype = serotype
        self.reference = references.get_reference(serotype)
        self.ref_seq = self.reference['sequence']

    def call_variants(self, query_sequence: str, ref_window: Optional[tuple] = None) -> List[Variant]:
        """
        Align query_sequence to the reference and call SNPs and indels.

        Args:
            query_sequence: Query genome sequence (need not be pre-aligned
                             or the same length as the reference).
            ref_window: Optional (start, end) 0-indexed slice of the
                        reference to align against instead of the full
                        genome; returned variant positions are still in
                        full-reference coordinates. Needed for short/partial
                        queries: plain global edit distance against the
                        *full* reference degenerates into scattering the
                        query across non-contiguous matching bases anywhere
                        in the genome once the length gap is large, because
                        with a 4-letter alphabet that's often cheaper than
                        aligning contiguously and eating a handful of real
                        mismatches. Callers analyzing fragments should
                        localize the query first (e.g. k-mer seeding) and
                        pass a padded window here; full-genome callers (the
                        default) are unaffected.

        Returns:
            List of Variant objects, in reference-coordinate order.
        """
        query_seq = query_sequence.upper()
        if ref_window:
            win_start, win_end = ref_window
            align_ref = self.ref_seq[win_start:win_end]
        else:
            win_start = 0
            align_ref = self.ref_seq

        ops = Levenshtein.editops(align_ref, query_seq)

        variants = []
        i, n = 0, len(ops)

        while i < n:
            op = ops[i]

            if op.tag == 'replace':
                ref_base = align_ref[op.src_pos]
                alt_base = query_seq[op.dest_pos]
                if ref_base != 'N' and alt_base != 'N':
                    variant = self._annotate_snp(win_start + op.src_pos, ref_base, alt_base)
                    if variant:
                        variants.append(variant)
                i += 1

            elif op.tag == 'delete':
                j = i + 1
                while (j < n and ops[j].tag == 'delete'
                       and ops[j].src_pos == ops[j - 1].src_pos + 1
                       and ops[j].dest_pos == ops[j - 1].dest_pos):
                    j += 1
                start, end = op.src_pos, ops[j - 1].src_pos + 1
                variants.append(self._annotate_indel(win_start + start, align_ref[start:end], ''))
                i = j

            elif op.tag == 'insert':
                j = i + 1
                while (j < n and ops[j].tag == 'insert'
                       and ops[j].src_pos == op.src_pos
                       and ops[j].dest_pos == ops[j - 1].dest_pos + 1):
                    j += 1
                alt = query_seq[op.dest_pos:ops[j - 1].dest_pos + 1]
                variants.append(self._annotate_indel(win_start + op.src_pos, '', alt))
                i = j

            else:
                i += 1

        return variants

    def _annotate_snp(self, position: int, ref_base: str, alt_base: str) -> Optional[Variant]:
        """Annotate a single SNP (0-indexed reference genome position)."""
        protein_name, protein_start = references.get_protein_at_position(self.serotype, position)

        if not protein_name:
            return Variant(
                position=position + 1,
                ref_allele=ref_base,
                alt_allele=alt_base,
                variant_type='SNP',
                impact='MODIFIER',
                protein=None,
                protein_position=None,
                codon_change=None,
                aa_change=None,
                is_synonymous=False,
                is_frameshift=False,
                is_stop_gained=False,
                is_stop_lost=False
            )

        codon_pos_in_protein = position - protein_start
        codon_index = codon_pos_in_protein // 3
        pos_in_codon = codon_pos_in_protein % 3

        codon_start = protein_start + codon_index * 3
        codon_end = codon_start + 3
        if codon_end > len(self.ref_seq):
            return None

        ref_codon = self.ref_seq[codon_start:codon_end]
        if len(ref_codon) < 3 or any(b not in 'ACGT' for b in ref_codon):
            return None
        alt_codon = ref_codon[:pos_in_codon] + alt_base + ref_codon[pos_in_codon + 1:]

        ref_aa = GENETIC_CODE.get(ref_codon, 'X')
        alt_aa = GENETIC_CODE.get(alt_codon, 'X')

        is_synonymous = (ref_aa == alt_aa)
        is_stop_gained = (alt_aa == '*' and ref_aa != '*')
        is_stop_lost = (ref_aa == '*' and alt_aa != '*')

        if is_stop_gained or is_stop_lost:
            impact = 'HIGH'
        elif ref_aa == 'X' or alt_aa == 'X':
            impact = 'MODERATE'
        elif is_synonymous:
            impact = 'LOW'
        else:
            impact = 'MODERATE'

        return Variant(
            position=position + 1,
            ref_allele=ref_base,
            alt_allele=alt_base,
            variant_type='SNP',
            impact=impact,
            protein=protein_name,
            protein_position=codon_index + 1,
            codon_change=f"{ref_codon}>{alt_codon}",
            aa_change=f"{ref_aa}{codon_index + 1}{alt_aa}",
            is_synonymous=is_synonymous,
            is_frameshift=False,
            is_stop_gained=is_stop_gained,
            is_stop_lost=is_stop_lost
        )

    def _annotate_indel(self, ref_position: int, ref_allele: str, alt_allele: str) -> Variant:
        """Annotate an insertion or deletion (0-indexed reference genome position)."""
        protein_name, protein_start = references.get_protein_at_position(self.serotype, ref_position)
        length_diff = len(alt_allele) - len(ref_allele)
        is_frameshift = bool(protein_name) and (length_diff % 3 != 0)

        if protein_name:
            protein_position = (ref_position - protein_start) // 3 + 1
            impact = 'HIGH' if is_frameshift else 'MODERATE'
        else:
            protein_position = None
            impact = 'MODIFIER'

        return Variant(
            position=ref_position + 1,
            ref_allele=ref_allele or '-',
            alt_allele=alt_allele or '-',
            variant_type='INDEL',
            impact=impact,
            protein=protein_name,
            protein_position=protein_position,
            codon_change=None,
            aa_change=None,
            is_synonymous=False,
            is_frameshift=is_frameshift,
            is_stop_gained=False,
            is_stop_lost=False
        )

    def annotate_with_epitopes(self, variant: Variant, epitope_data: List[dict]) -> dict:
        """
        Cross-reference a variant with known epitopes.

        Args:
            variant: Variant object
            epitope_data: List of epitope dicts with 'protein', 'start_pos', 'end_pos'
                          (protein-relative amino acid coordinates, 1-indexed)

        Returns:
            Enriched variant dict with an 'epitope_hits' list.
        """
        variant_dict = variant.to_dict()
        variant_dict['epitope_hits'] = []

        if not variant.protein or variant.protein_position is None:
            return variant_dict

        for epitope in epitope_data:
            if epitope.get('protein') != variant.protein:
                continue
            if epitope['start_pos'] <= variant.protein_position <= epitope['end_pos']:
                variant_dict['epitope_hits'].append({
                    'epitope_id': epitope.get('epitope_id'),
                    'epitope_type': epitope.get('epitope_type'),
                    'immune_response': epitope.get('immune_response'),
                })

        return variant_dict


if __name__ == "__main__":
    ref = references.get_reference(2)
    annotator = VariantAnnotator(serotype=2)
    # Self-test: reference vs. itself should yield zero variants.
    variants = annotator.call_variants(ref['sequence'])
    print(f"DENV-2 reference vs. itself: {len(variants)} variants (expected 0)")
