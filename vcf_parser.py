"""
VCF parser with tabix support, VEP/SnpEff annotation detection,
coordinate-based gene lookup, and variant prioritization.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Iterator

try:
    import cyvcf2
except ImportError:
    raise SystemExit("cyvcf2 is required: pip install cyvcf2")

from gene_panels import lookup_gene_by_position, GENE_COORDINATES
from gene_lookup import get_coordinates_batch, lookup_rsid_positions


@dataclass
class VariantRecord:
    gene: str
    chrom: str
    pos: int
    ref: str
    alt: str
    rsid: str
    zygosity: str           # heterozygous / homozygous_alt / hemizygous / unknown
    quality: float | None
    depth: int | None
    alt_allele_freq: float | None
    filter_status: str
    consequence: str        # SO term or "unknown"
    hgvs_c: str            # coding HGVS or ""
    hgvs_p: str            # protein HGVS or ""
    impact: str            # HIGH / MODERATE / LOW / MODIFIER / unknown
    clinvar_significance: str  # Pathogenic / Likely_pathogenic / VUS / etc. or ""
    gnomad_af: float | None
    panel: str = ""

    # Computed priority score (lower = more important)
    priority: int = field(default=0, init=False)

    def __post_init__(self) -> None:
        self.priority = self._compute_priority()

    def _compute_priority(self) -> int:
        sig = self.clinvar_significance.lower()
        if "pathogenic" in sig and "likely" not in sig and "conflicting" not in sig:
            return 0
        if "likely_pathogenic" in sig or "likely pathogenic" in sig:
            return 1
        if self.impact == "HIGH":
            return 2
        if self.impact == "MODERATE":
            af = self.gnomad_af if self.gnomad_af is not None else 1.0
            return 3 if af < 0.01 else 6
        if "uncertain" in sig or "vus" in sig:
            return 4
        af = self.gnomad_af if self.gnomad_af is not None else 1.0
        if af < 0.001:
            return 5
        return 7

    def to_dict(self) -> dict:
        return {
            "gene": self.gene,
            "variant": self._variant_label(),
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "rsid": self.rsid,
            "zygosity": self.zygosity,
            "quality": self.quality,
            "depth": self.depth,
            "alt_allele_freq": self.alt_allele_freq,
            "filter_status": self.filter_status,
            "consequence": self.consequence,
            "hgvs_c": self.hgvs_c,
            "hgvs_p": self.hgvs_p,
            "impact": self.impact,
            "clinvar_significance": self.clinvar_significance,
            "gnomad_af": self.gnomad_af,
        }

    def _variant_label(self) -> str:
        parts = []
        if self.hgvs_c:
            parts.append(self.hgvs_c)
        if self.hgvs_p:
            parts.append(f"({self.hgvs_p})")
        if not parts:
            parts.append(f"{self.ref}>{self.alt}")
        return " ".join(parts)


class VCFParser:
    MAX_VARIANTS_PER_PANEL = 100

    # Allele frequency cutoffs for inclusion without strong evidence
    AF_HIGH_IMPACT = 1.0    # always include HIGH impact
    AF_MODERATE = 0.01      # missense: include if gnomAD AF < 1%
    AF_LOW_IMPACT = 0.001   # other: include if gnomAD AF < 0.1%

    def __init__(self, vcf_path: str) -> None:
        if not os.path.exists(vcf_path):
            raise FileNotFoundError(f"VCF file not found: {vcf_path}")
        self.vcf_path = vcf_path
        self._vcf = cyvcf2.VCF(vcf_path)
        self.samples = self._vcf.samples
        self._annotation_format = self._detect_annotation_format()
        self._has_tabix = os.path.exists(vcf_path + ".tbi") or os.path.exists(
            vcf_path.replace(".vcf.gz", ".tbi")
        )
        self._vep_fields: list[str] = []
        if self._annotation_format == "vep":
            self._vep_fields = self._parse_vep_header()
        self.genome_build = self._detect_genome_build()
        # Chromosome names actually present in this VCF (from tabix/header)
        self._seqnames: set[str] = set(self._vcf.seqnames) if self._has_tabix else set()
        # "chr" if VCF uses chr-prefixed chroms, "" if bare, None if unknown
        self._chrom_prefix: str | None = self._detect_chrom_prefix()

    def _detect_annotation_format(self) -> str:
        raw = self._vcf.raw_header
        if "##INFO=<ID=CSQ," in raw:
            return "vep"
        if "##INFO=<ID=ANN," in raw:
            return "snpeff"
        return "none"

    def _parse_vep_header(self) -> list[str]:
        for line in self._vcf.raw_header.splitlines():
            if line.startswith("##INFO=<ID=CSQ,") and "Format:" in line:
                fmt_part = line.split("Format:")[1].strip().rstrip('">').strip()
                return [f.strip() for f in fmt_part.split("|")]
        return []

    def _detect_genome_build(self) -> str:
        """
        Detect whether the VCF is hg38/GRCh38 or hg19/GRCh37 from the header.
        Checks ##reference= lines first, then contig lengths for chr1, then
        any mention of build keywords in the raw header.
        Returns "hg38" or "hg19". Defaults to "hg38" if uncertain.
        """
        raw = self._vcf.raw_header

        # 1. Explicit ##reference= line
        for line in raw.splitlines():
            if not line.startswith("##reference="):
                continue
            ref_lower = line.lower()
            if "grch38" in ref_lower or "hg38" in ref_lower:
                return "hg38"
            if "grch37" in ref_lower or "hg19" in ref_lower or "b37" in ref_lower:
                return "hg19"
            # Some files have just "38" or "37" in the path
            if re.search(r"[_./]38[_./]", ref_lower) or ref_lower.endswith("38"):
                return "hg38"
            if re.search(r"[_./]37[_./]", ref_lower) or ref_lower.endswith("37"):
                return "hg19"

        # 2. Contig length for chr1 / 1
        #    GRCh38 chr1 = 248,956,422
        #    GRCh37 chr1 = 249,250,621
        for line in raw.splitlines():
            if not line.startswith("##contig="):
                continue
            if not re.search(r"ID=(chr)?1[,>]", line):
                continue
            m = re.search(r"length=(\d+)", line)
            if m:
                length = int(m.group(1))
                if length == 248956422:
                    return "hg38"
                if length == 249250621:
                    return "hg19"

        # 3. Any keyword in the full header
        if "GRCh38" in raw or "hg38" in raw:
            return "hg38"
        if "GRCh37" in raw or "hg19" in raw or "b37" in raw.lower():
            return "hg19"

        return "hg38"  # default assumption

    def get_diagnostics(self) -> dict:
        """Return diagnostic information about the VCF file for display in the UI."""
        prefix = self._chrom_prefix
        if prefix == "chr":
            chrom_display = "chr1 style (with prefix)"
        elif prefix == "":
            chrom_display = "1 style (no chr prefix)"
        else:
            chrom_display = "unknown (will try both)"

        # Sample first 3 variants to confirm the file is readable
        sample_variants: list[str] = []
        try:
            for v in cyvcf2.VCF(self.vcf_path):
                sample_variants.append(f"{v.CHROM}:{v.POS} {v.REF}>{v.ALT[0] if v.ALT else '.'}")
                if len(sample_variants) >= 3:
                    break
        except Exception:
            pass

        # Report a few actual seqnames so the user can see what's in the index
        seqnames_sample = sorted(self._seqnames)[:8] if self._seqnames else []

        return {
            "genome_build": self.genome_build,
            "chrom_prefix": prefix,
            "chrom_style": chrom_display,
            "annotation_format": self._annotation_format,
            "has_tabix": self._has_tabix,
            "num_samples": len(self.samples),
            "vep_fields": len(self._vep_fields),
            "seqnames_sample": seqnames_sample,
            "sample_variants": sample_variants,   # first 3 variants in the file
        }

    def _detect_chrom_prefix(self) -> str | None:
        """
        Determine whether this VCF uses 'chr'-prefixed chromosome names or bare names.
        Returns "chr", "", or None (unknown).
        Prefers seqnames (from tabix index) over header parsing.
        """
        if self._seqnames:
            if any(s.startswith("chr") for s in self._seqnames):
                return "chr"
            return ""
        # Fall back to header contig lines
        raw = self._vcf.raw_header
        if "##contig=<ID=chr" in raw:
            return "chr"
        if re.search(r"##contig=<ID=\d", raw):
            return ""
        return None  # cannot determine

    def _normalize_chrom(self, chrom: str) -> str:
        """Return the chromosome name normalised to 'chr' prefix (for internal lookups)."""
        if chrom.startswith("chr"):
            if chrom == "chrMT":
                return "chrM"
            return chrom
        # bare number or name
        if chrom in ("MT", "M"):
            return "chrM"
        return f"chr{chrom}"

    def _get_vcf_chrom(self, gene_chrom: str) -> list[str]:
        """
        Return ordered candidate chromosome names to try when querying the VCF.
        gene_chrom always has a 'chr' prefix (e.g. 'chr19', 'chrM').
        Puts the name matching the VCF's own naming style first.
        """
        # Derive bare form (strip 'chr')
        bare = gene_chrom[3:] if gene_chrom.startswith("chr") else gene_chrom
        # Handle MT aliases
        if gene_chrom in ("chrM", "chrMT"):
            if self._chrom_prefix == "chr":
                return ["chrM", "chrMT"]
            if self._chrom_prefix == "":
                return ["MT", "M"]
            return ["chrM", "MT", "chrMT", "M"]

        if self._chrom_prefix == "chr":
            return [gene_chrom]           # VCF definitely uses prefix
        if self._chrom_prefix == "":
            return [bare]                 # VCF definitely uses bare names
        # Unknown — try both, prefixed first
        return [gene_chrom, bare]

    def _parse_vep_annotation(self, variant: cyvcf2.Variant) -> dict:
        csq_raw = variant.INFO.get("CSQ")
        if not csq_raw or not self._vep_fields:
            return {}

        # CSQ may be a comma-separated list of transcripts; pick highest impact
        impact_order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
        best: dict = {}
        best_score = 99

        entries = csq_raw.split(",") if isinstance(csq_raw, str) else [str(csq_raw)]
        for entry in entries:
            vals = entry.split("|")
            ann = dict(zip(self._vep_fields, vals))
            score = impact_order.get(ann.get("IMPACT", ""), 99)
            if score < best_score:
                best_score = score
                best = ann

        gnomad_af: float | None = None
        for af_field in ("gnomADe_AF", "gnomAD_AF", "AF"):
            raw = best.get(af_field, "")
            if raw and raw not in (".", ""):
                try:
                    gnomad_af = float(raw)
                    break
                except ValueError:
                    pass

        clnsig = best.get("CLIN_SIG", "") or best.get("ClinVar_CLNSIG", "")

        return {
            "gene": best.get("SYMBOL", "") or best.get("Gene", ""),
            "consequence": best.get("Consequence", "unknown").split("&")[0],
            "hgvs_c": best.get("HGVSc", "").split(":")[-1],
            "hgvs_p": best.get("HGVSp", "").split(":")[-1].replace("%3D", "="),
            "impact": best.get("IMPACT", "unknown"),
            "clinvar_significance": clnsig.replace("&", "/"),
            "gnomad_af": gnomad_af,
        }

    def _parse_snpeff_annotation(self, variant: cyvcf2.Variant) -> dict:
        ann_raw = variant.INFO.get("ANN")
        if not ann_raw:
            return {}

        impact_order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
        best: dict = {}
        best_score = 99

        entries = ann_raw.split(",") if isinstance(ann_raw, str) else [str(ann_raw)]
        for entry in entries:
            fields = entry.split("|")
            if len(fields) < 10:
                continue
            impact = fields[2]
            score = impact_order.get(impact, 99)
            if score < best_score:
                best_score = score
                best = {
                    "gene": fields[3],
                    "consequence": fields[1].split("&")[0],
                    "hgvs_c": fields[9],
                    "hgvs_p": fields[10] if len(fields) > 10 else "",
                    "impact": impact,
                }

        clnsig = ""
        for key in ("CLNSIG", "ClinVar_CLNSIG", "CLIN_SIG"):
            val = variant.INFO.get(key)
            if val:
                clnsig = str(val)
                break

        gnomad_af: float | None = None
        for key in ("gnomAD_AF", "AF_popmax", "AF"):
            val = variant.INFO.get(key)
            if val is not None:
                try:
                    gnomad_af = float(val)
                    break
                except (ValueError, TypeError):
                    pass

        return {**best, "clinvar_significance": clnsig, "gnomad_af": gnomad_af}

    def _get_clinvar_info(self, variant: cyvcf2.Variant) -> str:
        for key in ("CLNSIG", "ClinVar_CLNSIG", "CLIN_SIG", "CLINSIG"):
            val = variant.INFO.get(key)
            if val:
                return str(val).replace("_", " ")
        return ""

    def _get_gnomad_af(self, variant: cyvcf2.Variant) -> float | None:
        for key in ("gnomAD_AF", "gnomADe_AF", "AF_popmax", "AF_gnomad", "AF"):
            val = variant.INFO.get(key)
            if val is not None:
                try:
                    f = float(val)
                    if 0.0 <= f <= 1.0:
                        return f
                except (ValueError, TypeError):
                    pass
        return None

    def _get_zygosity(self, variant: cyvcf2.Variant) -> str:
        if not self.samples:
            return "unknown"
        gt = variant.genotypes[0] if variant.genotypes else None
        if gt is None:
            return "unknown"
        alleles = [a for a in gt[:2] if a >= 0]
        if not alleles:
            return "unknown"
        ref_count = alleles.count(0)
        alt_count = len(alleles) - ref_count
        if alt_count == 0:
            return "homozygous_ref"
        if ref_count == 0 and alt_count >= 1:
            return "homozygous_alt" if len(alleles) == 2 else "hemizygous"
        return "heterozygous"

    def _get_allele_freq(self, variant: cyvcf2.Variant) -> float | None:
        if not self.samples or not variant.genotypes:
            return None
        gt = variant.genotypes[0]
        alleles = [a for a in gt[:2] if a >= 0]
        if not alleles:
            return None
        alt_count = sum(1 for a in alleles if a > 0)
        return alt_count / len(alleles)

    def _should_include(self, ann: dict, clnsig: str, gnomad_af: float | None) -> bool:
        """Return True if the variant passes clinical relevance filters."""
        sig = (clnsig or ann.get("clinvar_significance", "")).lower()
        if "pathogenic" in sig:
            return True
        impact = ann.get("impact", "unknown")
        if impact == "HIGH":
            return True
        af = gnomad_af if gnomad_af is not None else ann.get("gnomad_af")
        if isinstance(af, str):
            try:
                af = float(af)
            except ValueError:
                af = None
        if impact == "MODERATE":
            return af is None or af < self.AF_MODERATE
        if impact in ("LOW", "MODIFIER", "unknown"):
            return af is not None and af < self.AF_LOW_IMPACT
        return False

    def _build_record(
        self,
        variant: cyvcf2.Variant,
        gene: str,
        ann: dict,
        panel: str,
    ) -> VariantRecord:
        clnsig = ann.get("clinvar_significance") or self._get_clinvar_info(variant)
        gnomad_af = ann.get("gnomad_af") or self._get_gnomad_af(variant)
        chrom = variant.CHROM
        return VariantRecord(
            gene=gene or ann.get("gene", ""),
            chrom=chrom,
            pos=variant.POS,
            ref=variant.REF,
            alt=variant.ALT[0] if variant.ALT else ".",
            rsid=variant.ID or ".",
            zygosity=self._get_zygosity(variant),
            quality=variant.QUAL,
            depth=variant.INFO.get("DP"),
            alt_allele_freq=self._get_allele_freq(variant),
            filter_status="PASS" if variant.FILTER is None else str(variant.FILTER),
            consequence=ann.get("consequence", "unknown"),
            hgvs_c=ann.get("hgvs_c", ""),
            hgvs_p=ann.get("hgvs_p", ""),
            impact=ann.get("impact", "unknown"),
            clinvar_significance=clnsig,
            gnomad_af=gnomad_af,
            panel=panel,
        )

    def _annotate(self, variant: cyvcf2.Variant) -> dict:
        if self._annotation_format == "vep":
            return self._parse_vep_annotation(variant)
        if self._annotation_format == "snpeff":
            return self._parse_snpeff_annotation(variant)
        return {}

    def _iter_region(self, chrom: str, start: int, end: int) -> Iterator[cyvcf2.Variant]:
        """Yield variants in a chromosomal region using tabix if available."""
        if self._has_tabix:
            for candidate_chrom in self._get_vcf_chrom(chrom):
                try:
                    region = f"{candidate_chrom}:{start}-{end}"
                    yielded = False
                    for v in self._vcf(region):
                        yielded = True
                        yield v
                    # Only stop trying candidates when we actually got variants back.
                    # An empty result could mean "wrong chromosome name" (cyvcf2 does
                    # NOT raise an exception for unknown chrom names — it just returns
                    # an empty iterator), so we must not unconditionally return here.
                    if yielded:
                        return
                    # Zero variants: candidate name may be wrong — try the next one.
                except Exception:
                    continue
        else:
            # Full scan fallback — slow for WGS but works without index.
            # Open a fresh handle so we always scan from the start.
            norm = self._normalize_chrom(chrom)
            for v in cyvcf2.VCF(self.vcf_path):
                vc = self._normalize_chrom(v.CHROM)
                if vc == norm and start <= v.POS <= end:
                    yield v

    def filter_panel_variants(
        self,
        panel_name: str,
        gene_set: set[str],
    ) -> list[VariantRecord]:
        """
        Return up to MAX_VARIANTS_PER_PANEL variants for the given gene panel,
        prioritized by clinical significance.

        Two-pass strategy:
          Pass 1 — clinical filter (rare/high-impact variants)
          Pass 2 — curated rsID lookup (common variants with known significance)
                   The rsID list comes from PANEL_KEY_RSIDS in gene_panels.py.
                   This ensures common functional variants (e.g. ESR1 PvuII,
                   CYP2D6 *4, COMT Val158Met) are always reported for panels
                   where they are clinically relevant.
        """
        from gene_panels import PANEL_KEY_RSIDS

        records: list[VariantRecord] = []
        seen: set[tuple] = set()

        # ── Pass 1: clinical filter ──────────────────────────────────────────
        # Resolve coordinates for the entire gene set in one call,
        # using the correct build so hg19 VCFs get hg19 coords.
        coord_map = get_coordinates_batch(list(gene_set), build=self.genome_build)

        for gene in gene_set:
            coords = coord_map.get(gene)
            if coords is None:
                continue
            chrom, start, end = coords

            for variant in self._iter_region(chrom, start, end):
                if not variant.ALT:
                    continue

                ann = self._annotate(variant)

                # Determine gene name
                annotated_gene = ann.get("gene", "")
                if not annotated_gene:
                    annotated_gene = lookup_gene_by_position(
                        self._normalize_chrom(variant.CHROM), variant.POS
                    ) or gene

                clnsig = ann.get("clinvar_significance") or self._get_clinvar_info(variant)
                gnomad_af = ann.get("gnomad_af") or self._get_gnomad_af(variant)

                if not self._should_include(ann, clnsig, gnomad_af):
                    continue

                key = (variant.CHROM, variant.POS, variant.REF, tuple(variant.ALT))
                if key in seen:
                    continue
                seen.add(key)

                rec = self._build_record(variant, annotated_gene, ann, panel_name)
                records.append(rec)

        # ── Pass 2: curated rsID lookup ──────────────────────────────────────
        # For panels with known common-variant rsIDs (hormones, PGx, etc.),
        # look them up directly — bypasses frequency/impact filter entirely.
        key_rsids = PANEL_KEY_RSIDS.get(panel_name, [])
        if key_rsids:
            rsid_variants = self.get_variants_by_rsids(key_rsids)
            for rv in rsid_variants:
                k = (rv.chrom, rv.pos, rv.ref, rv.alt)
                if k in seen:
                    continue
                seen.add(k)
                rv.panel = panel_name
                records.append(rv)

        # Sort by priority (ClinVar P/LP first, then HIGH impact, etc.),
        # then by gene name for stable ordering within the same priority bucket.
        records.sort(key=lambda r: (r.priority, r.gene))
        return records[: self.MAX_VARIANTS_PER_PANEL]

    def get_variants_for_question(
        self,
        genes: list[str],
        rsids: list[str] | None = None,
    ) -> list[VariantRecord]:
        """
        Retrieve variants for a patient-specific question.
        1. Batch-resolve coordinates for ALL genes using the VCF's detected genome build.
        2. Use tabix region query for every gene that has coordinates.
        3. Fall back to annotation-based scan only for genes where coordinates cannot be found.
        """
        records: list[VariantRecord] = []
        seen: set[tuple] = set()

        # Resolve coordinates for all requested genes using the correct build
        coord_map = get_coordinates_batch(genes, build=self.genome_build)

        coord_genes = {g: coords for g, coords in coord_map.items() if coords is not None}
        no_coord_genes = [g for g, coords in coord_map.items() if coords is None]

        # Fast path: tabix region query for all genes with known coordinates
        for gene, (chrom, start, end) in coord_genes.items():
            for variant in self._iter_region(chrom, start, end):
                if not variant.ALT:
                    continue
                key = (variant.CHROM, variant.POS, variant.REF, tuple(variant.ALT))
                if key in seen:
                    continue
                seen.add(key)
                ann = self._annotate(variant)
                rec = self._build_record(variant, gene, ann, "query")
                records.append(rec)

        # Annotation-based fallback only for genes whose coordinates couldn't be found at all
        if no_coord_genes:
            annot_records = self.search_variants_by_gene_names(no_coord_genes, seen=seen)
            records.extend(annot_records)
            for r in annot_records:
                seen.add((r.chrom, r.pos, r.ref, r.alt))

        # If specific rsIDs are requested, do a targeted scan
        if rsids:
            rsid_set = set(rsids)
            found_rsids = {rec.rsid for rec in records}
            missing = rsid_set - found_rsids
            if missing:
                for variant in cyvcf2.VCF(self.vcf_path):
                    if variant.ID and variant.ID in missing:
                        ann = self._annotate(variant)
                        gene_name = ann.get("gene", "") or lookup_gene_by_position(
                            self._normalize_chrom(variant.CHROM), variant.POS
                        ) or "unknown"
                        key = (variant.CHROM, variant.POS, variant.REF, tuple(variant.ALT or []))
                        if key not in seen:
                            seen.add(key)
                            rec = self._build_record(variant, gene_name, ann, "query")
                            records.append(rec)

        records.sort(key=lambda r: (r.priority, r.gene))
        return records

    def search_variants_by_gene_names(
        self,
        gene_names: list[str],
        max_variants: int = 50,
        seen: set | None = None,
    ) -> list[VariantRecord]:
        """
        Scan the VCF for variants annotated to any of the given gene names.
        Works for VEP or SnpEff annotated VCFs. Returns empty list for unannotated VCFs.
        Applies the same clinical relevance filter as panel analysis.
        """
        if self._annotation_format == "none":
            return []

        gene_set = {g.upper() for g in gene_names}
        seen = seen or set()
        records: list[VariantRecord] = []

        # Always open a fresh reader so we scan from the start
        for variant in cyvcf2.VCF(self.vcf_path):
            if len(records) >= max_variants:
                break
            if not variant.ALT:
                continue

            ann = self._annotate(variant)
            annotated_gene = ann.get("gene", "").upper()
            if annotated_gene not in gene_set:
                continue

            clnsig = ann.get("clinvar_significance") or self._get_clinvar_info(variant)
            gnomad_af = ann.get("gnomad_af") or self._get_gnomad_af(variant)

            if not self._should_include(ann, clnsig, gnomad_af):
                continue

            key = (variant.CHROM, variant.POS, variant.REF, tuple(variant.ALT))
            if key in seen:
                continue
            seen.add(key)

            rec = self._build_record(variant, ann.get("gene", gene_names[0]), ann, "query")
            records.append(rec)

        records.sort(key=lambda r: (r.priority, r.gene))
        return records

    def get_variants_by_rsids(self, rsids: list[str]) -> list[VariantRecord]:
        """
        Look up specific variants by rsID, bypassing all clinical filters.
        Resolves each rsID to a GRCh38 position via Ensembl (cached), then uses
        tabix for fast retrieval. Falls back to full scan for unresolvable rsIDs.
        Returns the exact genotype regardless of allele frequency or impact.
        """
        rsids = [r.strip() for r in rsids if r.strip()]
        if not rsids:
            return []

        # Resolve rsID → (chrom, pos) via Ensembl + cache, using the correct build
        pos_map = lookup_rsid_positions(rsids, build=self.genome_build)
        records: list[VariantRecord] = []
        seen: set[tuple] = set()
        found_rsids: set[str] = set()

        for rsid, coord in pos_map.items():
            if coord is None:
                continue
            chrom, pos = coord
            # Query a tight 1-bp window around the position
            for variant in self._iter_region(chrom, pos, pos):
                if not variant.ALT:
                    continue
                # Match by position (and optionally rsID field)
                if variant.POS != pos:
                    continue
                key = (variant.CHROM, variant.POS, variant.REF, tuple(variant.ALT))
                if key in seen:
                    continue
                seen.add(key)
                ann = self._annotate(variant)
                gene = ann.get("gene", "") or lookup_gene_by_position(
                    self._normalize_chrom(variant.CHROM), variant.POS
                ) or rsid
                rec = self._build_record(variant, gene, ann, "query")
                # Attach the requested rsID if the VCF doesn't have it
                if rec.rsid == ".":
                    rec.rsid = rsid
                records.append(rec)
                found_rsids.add(rsid)

        # Fallback full scan for rsIDs that Ensembl couldn't resolve
        unresolved = [r for r in rsids if r not in pos_map or pos_map[r] is None]
        if unresolved:
            rsid_set = set(unresolved)
            for variant in cyvcf2.VCF(self.vcf_path):
                if not rsid_set:
                    break
                if variant.ID and variant.ID in rsid_set:
                    key = (variant.CHROM, variant.POS, variant.REF, tuple(variant.ALT or []))
                    if key not in seen:
                        seen.add(key)
                        ann = self._annotate(variant)
                        gene = ann.get("gene", "") or lookup_gene_by_position(
                            self._normalize_chrom(variant.CHROM), variant.POS
                        ) or variant.ID or "unknown"
                        rec = self._build_record(variant, gene, ann, "query")
                        records.append(rec)
                    rsid_set.discard(variant.ID)

        return records

    def get_unfiltered_variants(
        self,
        genes: list[str],
        max_per_gene: int = 300,
        pass_only: bool = True,
    ) -> list[VariantRecord]:
        """
        Return ALL variants in gene regions with NO clinical significance filter.
        Includes common SNPs (any MAF), synonymous, intron variants — everything.
        Use for polygenic/GWAS analysis, pharmacogenomics star-alleles, etc.
        Results are capped at max_per_gene per gene to avoid flooding the response.
        """
        coord_map = get_coordinates_batch(genes, build=self.genome_build)
        records: list[VariantRecord] = []
        seen: set[tuple] = set()

        for gene, coords in coord_map.items():
            if coords is None:
                continue
            chrom, start, end = coords
            count = 0
            for variant in self._iter_region(chrom, start, end):
                if count >= max_per_gene:
                    break
                if not variant.ALT:
                    continue
                if pass_only and variant.FILTER is not None:
                    continue
                key = (variant.CHROM, variant.POS, variant.REF, tuple(variant.ALT))
                if key in seen:
                    continue
                seen.add(key)
                ann = self._annotate(variant)
                rec = self._build_record(variant, gene, ann, "query")
                records.append(rec)
                count += 1

        records.sort(key=lambda r: (r.gene, r.pos))
        return records

    @property
    def has_annotations(self) -> bool:
        """True if the VCF has VEP or SnpEff annotations."""
        return self._annotation_format != "none"
