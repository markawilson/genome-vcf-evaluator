"""
Real-time clinical database lookups for genomic variant analysis.

Overcomes LLM training-data recall bias by querying live sources:
  - NCBI ClinVar  — Pathogenic / Likely Pathogenic variants per gene
  - EBI GWAS Catalog — genome-wide significant associations per gene
  - PharmGKB       — pharmacogenomic clinical annotations per gene

Results are cached on disk (24 h TTL) so repeated queries within a session
are near-instant without hammering external APIs.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Optional

import requests

# ── Configuration ─────────────────────────────────────────────────────────────

CACHE_DIR = Path.home() / ".cache" / "genome-vcf-evaluator" / "db_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
CACHE_TTL_HOURS = 24.0   # refresh after 24 hours

NCBI_ESEARCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
GWAS_API_BASE    = "https://www.ebi.ac.uk/gwas/rest/api"
PHARMGKB_API_BASE = "https://api.pharmgkb.org/v1"

_HEADERS = {
    "Accept": "application/json",
    "User-Agent": "genome-vcf-evaluator/1.0 (research; not-for-clinical-use)",
}

# ── Disk-cache helpers ─────────────────────────────────────────────────────────

def _cache_path(prefix: str, key: str) -> Path:
    h = hashlib.md5(key.encode()).hexdigest()[:10]
    return CACHE_DIR / f"{prefix}_{h}.json"


def _read_cache(path: Path) -> Optional[list]:
    if not path.exists():
        return None
    age = (time.time() - path.stat().st_mtime) / 3600.0
    if age > CACHE_TTL_HOURS:
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _write_cache(path: Path, data: list) -> None:
    try:
        path.write_text(json.dumps(data, default=str))
    except Exception:
        pass


# ── ClinVar (NCBI E-utilities) ─────────────────────────────────────────────────

def query_clinvar(gene: str, max_results: int = 25) -> list[dict]:
    """
    Return Pathogenic / Likely Pathogenic ClinVar variants for *gene*.

    Each item:
      rsid, clinvar_id, name (HGVS title), clinical_significance,
      conditions (list[str]), gene, source="ClinVar"

    Set env var NCBI_API_KEY to raise rate limit from 3 → 10 req/s.
    """
    gene = gene.upper()
    cache = _cache_path("clinvar", gene)
    cached = _read_cache(cache)
    if cached is not None:
        return cached

    api_key = os.environ.get("NCBI_API_KEY", "")
    delay = 0.12 if api_key else 0.4   # NCBI rate-limit safe interval

    try:
        # ── Step 1: search ClinVar for P / LP variant IDs ──
        params: dict = {
            "db": "clinvar",
            "term": (
                f"{gene}[gene] AND "
                '("pathogenic"[clinsig] OR "likely pathogenic"[clinsig])'
            ),
            "retmax": max_results,
            "retmode": "json",
        }
        if api_key:
            params["api_key"] = api_key

        r = requests.get(NCBI_ESEARCH, params=params, headers=_HEADERS, timeout=12)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])

        if not ids:
            _write_cache(cache, [])
            return []

        time.sleep(delay)

        # ── Step 2: fetch summaries ──
        params2: dict = {
            "db": "clinvar",
            "id": ",".join(ids),
            "retmode": "json",
        }
        if api_key:
            params2["api_key"] = api_key

        r2 = requests.get(NCBI_ESUMMARY, params=params2, headers=_HEADERS, timeout=15)
        r2.raise_for_status()
        result = r2.json().get("result", {})

        variants: list[dict] = []
        for vid in ids:
            entry = result.get(vid)
            if not entry or entry.get("error"):
                continue

            # Try to pull rsID from variation_set → variation → xrefs
            rsid = None
            for vset in entry.get("variation_set", []):
                for var in vset.get("variation", []):
                    for xref in var.get("xrefs", []):
                        if xref.get("db", "").lower() in ("dbsnp", "rs"):
                            rsid = f"rs{xref['id']}"
                            break

            sig = ""
            cs = entry.get("clinical_significance")
            if isinstance(cs, dict):
                sig = cs.get("description", "")
            elif isinstance(cs, str):
                sig = cs

            conditions = [
                t.get("trait_name", "")
                for t in entry.get("trait_set", [])
                if t.get("trait_name")
            ]

            variants.append({
                "rsid": rsid,
                "clinvar_id": vid,
                "name": entry.get("title", ""),
                "clinical_significance": sig,
                "conditions": conditions[:4],
                "gene": gene,
                "source": "ClinVar",
            })

        _write_cache(cache, variants)
        return variants

    except Exception as exc:
        return [{"error": str(exc), "source": "ClinVar", "gene": gene}]


# ── GWAS Catalog (EBI) ─────────────────────────────────────────────────────────

# Functional-class priority (lower = more likely to be directly causal / interpretable)
_FC_PRIORITY: dict = {
    "stop_gained": 0, "stop_lost": 0, "start_lost": 0,
    "frameshift_variant": 1, "splice_acceptor_variant": 1, "splice_donor_variant": 1,
    "missense_variant": 2, "inframe_insertion": 3, "inframe_deletion": 3,
    "synonymous_variant": 4, "splice_region_variant": 5,
}

def query_gwas_catalog(gene: str, max_rsids: int = 60) -> list[dict]:
    """
    Return GWAS-associated SNPs for *gene* from EBI GWAS Catalog.

    Strategy: one fast call to the `findByGene` SNP endpoint (1-2 s), which
    returns all SNPs in the GWAS Catalog that map to the gene.  We de-duplicate
    by rsID, sort coding/splice variants first, and return the rsID list.

    Claude uses these rsIDs to query the patient VCF in rsID mode.
    Claude's own training knowledge interprets whatever is found.

    Each item:
      rsid, functional_class, gene, source="GWAS Catalog"
    A trailing summary item is appended when the catalog contains more SNPs
    than we return.
    """
    gene = gene.upper()
    cache = _cache_path("gwas", gene)
    cached = _read_cache(cache)
    if cached is not None:
        return cached

    try:
        url = f"{GWAS_API_BASE}/singleNucleotidePolymorphisms/search/findByGene"
        # Fetch up to 100 at once (max page size the API supports).
        # For genes with many GWAS hits we'll note how many total there are.
        params = {"geneName": gene, "size": 100}
        r = requests.get(url, params=params, headers=_HEADERS, timeout=15)
        r.raise_for_status()

        data = r.json()
        snps = data.get("_embedded", {}).get("singleNucleotidePolymorphisms", [])
        total = data.get("page", {}).get("totalElements", len(snps))

        seen: set = set()
        results: list[dict] = []

        for snp in snps:
            rs = snp.get("rsId")
            if not rs or rs in seen:
                continue
            seen.add(rs)

            fc = snp.get("functionalClass", "intergenic_variant")
            prio = _FC_PRIORITY.get(fc, 10)

            results.append({
                "rsid": rs,
                "functional_class": fc,
                "_prio": prio,
                "gene": gene,
                "source": "GWAS Catalog",
            })

        # Coding/splice variants first, then others in original order
        results.sort(key=lambda x: x["_prio"])
        for item in results:
            del item["_prio"]

        results = results[:max_rsids]

        # Append a summary note when the catalog has more than we returned
        n_unique = len(seen)
        if total > n_unique or n_unique > max_rsids:
            results.append({
                "summary": (
                    f"GWAS Catalog contains {total} associations "
                    f"({n_unique} unique rsIDs) for {gene}. "
                    f"Returning top {len(results)} prioritised by functional impact. "
                    "Use query_vcf in rsID mode with the listed rsIDs to check the patient's genotypes."
                ),
                "source": "GWAS Catalog",
                "gene": gene,
            })

        _write_cache(cache, results)
        return results

    except Exception as exc:
        return [{"error": str(exc), "source": "GWAS Catalog", "gene": gene}]


# ── PharmGKB ──────────────────────────────────────────────────────────────────

import re as _re

# Evidence level ordering: 1A is strongest, 4 is weakest
_PGKB_LEVEL_ORDER = {"1A": 0, "1B": 1, "2A": 2, "2B": 3, "3": 4, "4": 5}

def query_pharmgkb(gene: str, max_results: int = 25) -> list[dict]:
    """
    Return pharmacogenomic clinical annotations from PharmGKB for *gene*.

    Uses the `location.genes.symbol` filter (not `gene.symbol`).
    Results sorted by evidence level (1A strongest → 4 weakest).

    Each item:
      rsids (list), haplotypes (list), drug, phenotype, evidence_level,
      annotation_name, gene, source="PharmGKB"
    """
    gene = gene.upper()
    cache = _cache_path("pharmgkb", gene)
    cached = _read_cache(cache)
    if cached is not None:
        return cached

    try:
        url = f"{PHARMGKB_API_BASE}/data/clinicalAnnotation"
        params = {"location.genes.symbol": gene}
        r = requests.get(url, params=params, headers=_HEADERS, timeout=15)
        r.raise_for_status()

        annotations: list[dict] = []
        for item in r.json().get("data", []):
            name = item.get("name", "")

            # Extract rsIDs and star-allele haplotypes from the annotation name
            rsids = _re.findall(r"rs\d+", name)
            haplotypes = _re.findall(r"\*\d+\w*", name)

            drugs = [c.get("name", "") for c in item.get("relatedChemicals", [])]
            drug_str = ", ".join(d for d in drugs if d) or "—"

            phenos = [d.get("name", "") for d in item.get("relatedDiseases", [])]
            pheno_str = "; ".join(p for p in phenos if p)

            level_obj = item.get("levelOfEvidence", {})
            level = (
                level_obj.get("term", "")
                if isinstance(level_obj, dict)
                else str(level_obj)
            )

            annotations.append({
                "rsids": rsids,
                "haplotypes": haplotypes[:6],  # cap to avoid oversized payloads
                "drug": drug_str,
                "phenotype": pheno_str,
                "evidence_level": level,
                "annotation_name": name[:120],
                "_sort": _PGKB_LEVEL_ORDER.get(level, 99),
                "gene": gene,
                "source": "PharmGKB",
            })

        # Sort by evidence strength, then trim
        annotations.sort(key=lambda x: x["_sort"])
        for a in annotations:
            del a["_sort"]
        annotations = annotations[:max_results]

        _write_cache(cache, annotations)
        return annotations

    except Exception as exc:
        return [{"error": str(exc), "source": "PharmGKB", "gene": gene}]


# ── Combined entry point ───────────────────────────────────────────────────────

def lookup_gene_databases(
    genes: list[str],
    databases: Optional[list[str]] = None,
) -> dict:
    """
    Run all requested database queries for up to 5 genes.

    databases: any subset of ["clinvar", "gwas", "pharmgkb"]
               (default: all three)

    Returns:
        {
          "GENE1": {"clinvar": [...], "gwas_catalog": [...], "pharmgkb": [...]},
          "GENE2": {...},
          ...
        }
    """
    if databases is None:
        databases = ["clinvar", "gwas", "pharmgkb"]

    results: dict = {}
    for gene in genes[:5]:
        g = gene.upper()
        gr: dict = {}
        if "clinvar" in databases:
            gr["clinvar"] = query_clinvar(g)
            time.sleep(0.35)          # NCBI rate limit
        if "gwas" in databases:
            gr["gwas_catalog"] = query_gwas_catalog(g)
        if "pharmgkb" in databases:
            gr["pharmgkb"] = query_pharmgkb(g)
        results[g] = gr

    return results
