"""
Dynamic hg38 gene coordinate lookup with local disk cache.

Priority order for every gene symbol:
  1. Hardcoded GENE_COORDINATES (instant, offline)
  2. Local disk cache  ~/.genome_vcf_evaluator/gene_coords_cache.json  (instant, offline)
  3. mygene.info REST API  (requires internet, ~200 ms per batch)

Results from step 3 are written to the cache so subsequent runs are instant.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Optional

from gene_panels import GENE_COORDINATES

CACHE_DIR = Path.home() / ".genome_vcf_evaluator"
CACHE_FILE = CACHE_DIR / "gene_coords_cache.json"           # hg38
CACHE_FILE_HG19 = CACHE_DIR / "gene_coords_cache_hg19.json"
RSID_CACHE_FILE = CACHE_DIR / "rsid_coords_cache.json"      # GRCh38
RSID_CACHE_FILE_HG19 = CACHE_DIR / "rsid_coords_cache_hg19.json"
MYGENE_URL = "https://mygene.info/v3/query"
ENSEMBL_VARIATION_URL = "https://rest.ensembl.org/variation/homo_sapiens"
ENSEMBL_VARIATION_URL_HG19 = "https://grch37.rest.ensembl.org/variation/homo_sapiens"
REQUEST_TIMEOUT = 10  # seconds


# ── Cache helpers ──────────────────────────────────────────────────────────────

def _load_cache() -> dict:
    if CACHE_FILE.exists():
        try:
            return json.loads(CACHE_FILE.read_text(encoding="utf-8"))
        except Exception:
            pass
    return {}


def _save_cache(cache: dict) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_FILE.write_text(json.dumps(cache, indent=2), encoding="utf-8")


# ── mygene.info queries ────────────────────────────────────────────────────────

def _parse_genomic_pos(pos) -> Optional[tuple[str, int, int]]:
    """Parse a mygene genomic_pos entry (dict or first element of list)."""
    if isinstance(pos, list):
        pos = pos[0] if pos else None
    if not isinstance(pos, dict):
        return None
    chrom = str(pos.get("chr", "")).strip()
    start = pos.get("start")
    end = pos.get("end")
    if not chrom or start is None or end is None:
        return None
    try:
        start, end = int(start), int(end)
    except (TypeError, ValueError):
        return None
    if start <= 0 or end <= 0:
        return None
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    if chrom == "chrMT":
        chrom = "chrM"
    return (chrom, start, end)


def _query_mygene_batch(
    genes: list[str], build: str = "hg38"
) -> dict[str, Optional[tuple[str, int, int]]]:
    """
    Batch-query mygene.info for gene coordinates.
    build: "hg38" → genomic_pos field; "hg19" → genomic_pos_hg19 field.
    Returns {gene_symbol: (chrom, start, end) or None}.
    """
    if not genes:
        return {}

    pos_field = "genomic_pos_hg19" if build == "hg19" else "genomic_pos"
    results: dict[str, Optional[tuple[str, int, int]]] = {g: None for g in genes}

    try:
        body = urllib.parse.urlencode(
            {
                "q": ",".join(genes),
                "scopes": "symbol",
                "species": "human",
                "fields": f"{pos_field},symbol",
                "size": str(min(len(genes) * 3, 1000)),
            }
        ).encode()

        req = urllib.request.Request(
            MYGENE_URL,
            data=body,
            headers={
                "User-Agent": "GenomeVCFAnalyzer/1.0",
                "Content-Type": "application/x-www-form-urlencoded",
            },
        )
        with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as resp:
            hits = json.loads(resp.read())

        if isinstance(hits, dict):
            hits = hits.get("hits", [])

        for hit in hits:
            symbol = str(hit.get("symbol") or hit.get("query") or "").upper().strip()
            if symbol not in results:
                continue
            coords = _parse_genomic_pos(hit.get(pos_field))
            if coords and results[symbol] is None:
                results[symbol] = coords

    except (urllib.error.URLError, OSError, json.JSONDecodeError):
        pass  # offline or API error — return None for all

    return results


# ── Public API ─────────────────────────────────────────────────────────────────

def get_coordinates(gene: str, build: str = "hg38") -> Optional[tuple[str, int, int]]:
    """Return (chrom, start, end) for a gene in the requested build, or None."""
    gene = gene.upper().strip()

    # 1. Hardcoded database (hg38 only)
    if build == "hg38" and gene in GENE_COORDINATES:
        return GENE_COORDINATES[gene]

    # 2. Disk cache
    cache_file = CACHE_FILE_HG19 if build == "hg19" else CACHE_FILE
    cache = json.loads(cache_file.read_text()) if cache_file.exists() else {}
    if gene in cache:
        entry = cache[gene]
        return tuple(entry) if entry else None  # type: ignore[return-value]

    # 3. Online lookup (single gene)
    result = _query_mygene_batch([gene], build=build)
    coords = result.get(gene)
    cache[gene] = list(coords) if coords else None
    _save_cache(cache)
    return coords


def get_coordinates_batch(
    genes: list[str],
    build: str = "hg38",
) -> dict[str, Optional[tuple[str, int, int]]]:
    """
    Return coordinates for a list of genes in the requested genome build.
    Priority: hardcoded DB (hg38 only) → build-specific disk cache → mygene.info.
    Returns {GENE: (chrom, start, end) or None}.
    """
    genes = [g.upper().strip() for g in genes]
    cache_file = CACHE_FILE_HG19 if build == "hg19" else CACHE_FILE
    cache: dict = {}
    if cache_file.exists():
        try:
            cache = json.loads(cache_file.read_text(encoding="utf-8"))
        except Exception:
            pass

    results: dict[str, Optional[tuple[str, int, int]]] = {}
    need_online: list[str] = []

    for gene in genes:
        # Hardcoded database covers hg38 only
        if build == "hg38" and gene in GENE_COORDINATES:
            results[gene] = GENE_COORDINATES[gene]
        elif gene in cache:
            entry = cache[gene]
            results[gene] = tuple(entry) if entry else None  # type: ignore[assignment]
        else:
            need_online.append(gene)

    if need_online:
        online = _query_mygene_batch(need_online, build=build)
        for gene, coords in online.items():
            results[gene] = coords
            cache[gene] = list(coords) if coords else None
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(json.dumps(cache, indent=2), encoding="utf-8")

    return results


def cache_stats() -> dict:
    """Return info about the local cache for display."""
    cache = _load_cache()
    found = sum(1 for v in cache.values() if v is not None)
    return {
        "total_cached": len(cache),
        "with_coordinates": found,
        "not_found": len(cache) - found,
        "cache_path": str(CACHE_FILE),
    }


# ── rsID → genomic position lookup ────────────────────────────────────────────

def _load_rsid_cache() -> dict:
    if RSID_CACHE_FILE.exists():
        try:
            return json.loads(RSID_CACHE_FILE.read_text(encoding="utf-8"))
        except Exception:
            pass
    return {}


def _save_rsid_cache(cache: dict) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    RSID_CACHE_FILE.write_text(json.dumps(cache, indent=2), encoding="utf-8")


def _query_ensembl_rsids(
    rsids: list[str], build: str = "hg38"
) -> dict[str, Optional[tuple[str, int]]]:
    """
    Batch-query Ensembl REST API for positions of rsIDs.
    build: "hg38" → GRCh38 assembly via rest.ensembl.org
           "hg19" → GRCh37 assembly via grch37.rest.ensembl.org
    Returns {rsid: (chrom, pos)} or {rsid: None} for not found.
    Ensembl accepts up to 200 IDs per POST.
    """
    results: dict[str, Optional[tuple[str, int]]] = {r: None for r in rsids}
    assembly_name = "GRCh37" if build == "hg19" else "GRCh38"
    url = ENSEMBL_VARIATION_URL_HG19 if build == "hg19" else ENSEMBL_VARIATION_URL

    # Process in batches of 200
    for i in range(0, len(rsids), 200):
        batch = rsids[i : i + 200]
        try:
            body = json.dumps({"ids": batch}).encode()
            req = urllib.request.Request(
                url,
                data=body,
                headers={
                    "Content-Type": "application/json",
                    "Accept": "application/json",
                    "User-Agent": "GenomeVCFAnalyzer/1.0",
                },
            )
            with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as resp:
                data = json.loads(resp.read())

            for rsid, info in data.items():
                if not isinstance(info, dict):
                    continue
                for mapping in info.get("mappings", []):
                    if mapping.get("assembly_name") != assembly_name:
                        continue
                    loc = mapping.get("location", "")   # e.g. "22:19951271-19951271"
                    parts = loc.split(":")
                    if len(parts) != 2:
                        continue
                    raw_chrom = parts[0]
                    pos_str = parts[1].split("-")[0]
                    try:
                        pos = int(pos_str)
                    except ValueError:
                        continue
                    chrom = raw_chrom if raw_chrom.startswith("chr") else f"chr{raw_chrom}"
                    if chrom == "chrMT":
                        chrom = "chrM"
                    results[rsid] = (chrom, pos)
                    break  # take first matching assembly mapping

        except (urllib.error.URLError, OSError, json.JSONDecodeError):
            pass  # offline or API error

    return results


def lookup_rsid_positions(
    rsids: list[str], build: str = "hg38"
) -> dict[str, Optional[tuple[str, int]]]:
    """
    Return (chrom, pos) for each rsID in the requested genome build.
    Checks build-specific disk cache first; queries Ensembl for uncached rsIDs.
    Returns {rsid: (chrom, pos)} or {rsid: None}.
    """
    cache_file = RSID_CACHE_FILE_HG19 if build == "hg19" else RSID_CACHE_FILE
    cache: dict = {}
    if cache_file.exists():
        try:
            cache = json.loads(cache_file.read_text(encoding="utf-8"))
        except Exception:
            pass

    results: dict[str, Optional[tuple[str, int]]] = {}
    need_online: list[str] = []

    for rsid in rsids:
        if rsid in cache:
            entry = cache[rsid]
            results[rsid] = tuple(entry) if entry else None  # type: ignore[assignment]
        else:
            need_online.append(rsid)

    if need_online:
        online = _query_ensembl_rsids(need_online, build=build)
        for rsid, pos_tuple in online.items():
            results[rsid] = pos_tuple
            cache[rsid] = list(pos_tuple) if pos_tuple else None
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(json.dumps(cache, indent=2), encoding="utf-8")

    return results
