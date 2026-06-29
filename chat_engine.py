"""
Agentic chat engine — Claude with two tools it can call autonomously:
  • query_databases — live ClinVar / GWAS Catalog / PharmGKB lookup
  • query_vcf       — patient VCF genotype query

Flow per turn:
  1. Send user message + tool definitions to Claude
  2. Claude may call query_databases to discover relevant rsIDs from live sources
  3. Claude may call query_vcf (returns variants from the local VCF)
  4. Steps 2-3 can interleave; Claude stops when it has enough to answer
  5. Final text response is yielded as a stream of events for the UI
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Generator, List, Union

import anthropic

from database_lookup import lookup_gene_databases
from gene_lookup import get_coordinates_batch
from vcf_parser import VCFParser, VariantRecord

# ── Tool definition ────────────────────────────────────────────────────────────

QUERY_VCF_TOOL: dict = {
    "name": "query_vcf",
    "description": (
        "Query the patient's whole-genome VCF file for variants. Supports three modes:\n\n"
        "1. CLINICAL MODE (default): genes + clinical filter — returns only ClinVar P/LP, "
        "HIGH-impact, rare missense, and ultra-rare variants. Best for Mendelian disease.\n\n"
        "2. UNFILTERED MODE (include_common=true): genes + NO filter — returns ALL variants "
        "in the gene regions including common SNPs (any MAF). Use for polygenic traits, "
        "GWAS loci, pharmacogenomics star-alleles, hormonal variants (ESR1 PvuII, etc.), "
        "or any time you need genotypes at common positions.\n\n"
        "3. RSID MODE: rsids list — looks up specific variants by rsID (e.g. rs4680 for "
        "COMT Val158Met, rs9340799 for ESR1 PvuII) regardless of frequency or impact. "
        "Use when you know the exact rsID of a clinically relevant SNP."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "genes": {
                "type": "array",
                "items": {"type": "string"},
                "description": (
                    "Gene symbols to query, uppercase (e.g. ['ESR1', 'ESR2', 'PGR']). "
                    "Required for clinical and unfiltered modes. Optional if rsids provided."
                ),
            },
            "rsids": {
                "type": "array",
                "items": {"type": "string"},
                "description": (
                    "Specific rsIDs to look up (e.g. ['rs4680', 'rs9340799', 'rs2234693']). "
                    "Returns exact genotype at that position, bypassing all frequency and impact "
                    "filters. Use for known GWAS/pharmacogenomic SNPs where you know the rsID."
                ),
            },
            "include_common": {
                "type": "boolean",
                "description": (
                    "If true, bypass ALL frequency and impact filters and return every variant "
                    "in the gene regions (capped at 300 per gene). Use for polygenic risk, "
                    "hormonal SNPs, GWAS loci, or any common variant analysis. "
                    "Defaults to false (clinical filter only)."
                ),
            },
            "reason": {
                "type": "string",
                "description": "One sentence explaining what you are looking for and why.",
            },
        },
        "required": ["reason"],
    },
}

QUERY_DATABASES_TOOL: dict = {
    "name": "query_databases",
    "description": (
        "Query live clinical databases to discover which variants are known to be "
        "relevant for a gene — BEFORE or alongside querying the patient VCF.\n\n"
        "This overcomes training-data recall bias: instead of relying only on what you "
        "remember about a gene's variants, you fetch the current published record.\n\n"
        "Three sources (specify any subset via 'databases', or omit for all three):\n"
        "• 'clinvar'  — NCBI ClinVar: all Pathogenic/Likely Pathogenic variants for the gene\n"
        "• 'gwas'     — EBI GWAS Catalog: genome-wide significant associations (p≤10⁻⁶)\n"
        "               Use for: hormonal traits, mood, cognition, metabolism, pain, migraine\n"
        "• 'pharmgkb' — PharmGKB: drug-variant-phenotype annotations\n"
        "               Use for: antidepressant response, PGx, drug metabolism\n\n"
        "After calling this tool, use the returned rsIDs to call query_vcf in rsID mode "
        "to look up the patient's actual genotype at those positions."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "genes": {
                "type": "array",
                "items": {"type": "string"},
                "description": (
                    "Gene symbols to query (uppercase, max 5 per call). "
                    "For broad topics, call this tool once per gene or small batch."
                ),
            },
            "databases": {
                "type": "array",
                "items": {"type": "string", "enum": ["clinvar", "gwas", "pharmgkb"]},
                "description": (
                    "Which databases to query. Omit to query all three. "
                    "Recommended subsets: ['gwas'] for GWAS/polygenic/hormonal/mood traits; "
                    "['clinvar'] for Mendelian disease genes; "
                    "['pharmgkb'] for drug-metabolism genes (CYP2D6, CYP2C19, etc.)."
                ),
            },
            "reason": {
                "type": "string",
                "description": "One sentence explaining what you are looking for and why.",
            },
        },
        "required": ["genes", "reason"],
    },
}

# ── System prompt ──────────────────────────────────────────────────────────────

CHAT_SYSTEM_PROMPT = """\
You are an expert clinical genomicist with real-time access to:
  1. Live clinical databases (query_databases tool)
  2. The patient's whole-genome VCF (query_vcf tool)

## RECOMMENDED WORKFLOW

For any non-trivial question, use this two-step approach:

**Step 1 — query_databases first** (unless you are 100% confident you already know \
all relevant rsIDs for the question):
  - Fetches live ClinVar P/LP variants, GWAS associations, and/or PharmGKB PGx data
  - Returns the current published list of clinically important rsIDs for those genes
  - Overcomes the critical limitation that your training memory only knows \
"famous" variants and misses newer or less-cited GWAS hits
  - Use databases=["gwas"] for hormonal/mood/polygenic traits
  - Use databases=["clinvar"] for Mendelian disease genes
  - Use databases=["pharmgkb"] for drug metabolism genes
  - Omit databases to query all three

**Step 2 — query_vcf with the discovered rsIDs**:
  - Take the rsIDs returned from Step 1 and look them up in rsID mode
  - Also look up any additional rsIDs you already know are important
  - This gives you a comprehensive genotype picture grounded in current evidence

Skip Step 1 only when:
  - You are asked about a specific rsID you already have
  - You are doing ancestry-informative marker lookups (rsIDs listed below)
  - The question is purely clinical/Mendelian and you are very confident in the rsID list

## query_vcf — three modes

**Clinical mode** (default): rare/high-impact variants only. Best for Mendelian disease.

**Unfiltered mode** (include_common=true): ALL variants including common SNPs. Use for:
  - Hormonal/receptor variants, polygenic GWAS loci, neurotransmitter SNPs

**rsID mode** (rsids=[...]): Direct lookup of specific SNPs. Use after query_databases \
returns rsIDs, or for well-known SNPs you already have memorized:
  - COMT Val158Met → rs4680
  - MTHFR C677T → rs1801133, A1298C → rs1801131
  - Factor V Leiden → rs6025, Prothrombin G20210A → rs1799963
  - ESR1 PvuII → rs2234693, XbaI → rs9340799

**Ancestry markers** (rsid mode):
  - SLC24A5 (European skin) → rs1426654 | SLC45A2 → rs16891982
  - HERC2/OCA2 (blue eyes) → rs12913832 | MC1R (red hair) → rs1805007/8/9
  - EDAR (East Asian) → rs3827760 | LCT (lactase) → rs4988235
  - ALDH2 (alcohol flush) → rs671 | ADH1B → rs1229984
  - ACKR1 Duffy null (African) → rs2814778
  - APOE ε2/ε4 → rs7412, rs429358 | APOL1 G1 → rs73885319, rs60910145

## Interpretation principles
- Integrate zygosity, population frequency, effect size (OR/β from GWAS), and ClinVar \
  classification into every interpretation
- For GWAS/polygenic findings: note that individual SNPs have modest effects (OR 1.1–1.3) \
  and population-level risk does not directly translate to individual prediction
- Always report what the reference genotype means when no variant is found

IMPORTANT DISCLAIMER: All findings are for research and educational purposes only. \
They do not constitute medical advice. Clinical decisions must be made by qualified \
healthcare professionals using CLIA-certified laboratory testing.\
"""

# ── Event types yielded to the UI ─────────────────────────────────────────────

@dataclass
class DatabaseLookupEvent:
    """Fired when Claude calls query_databases (before the results arrive)."""
    genes: List[str]
    databases: List[str]
    reason: str

@dataclass
class DatabaseResultEvent:
    """Fired when query_databases results are ready."""
    genes: List[str]
    databases: List[str]
    counts: dict    # {gene: {source: n_results}}
    rsids_found: List[str]   # all rsIDs discovered (for UI summary)

@dataclass
class ToolCallEvent:
    genes: List[str]
    reason: str
    fetched_online: List[str]     # genes whose coords were looked up via mygene.info

@dataclass
class ToolResultEvent:
    genes: List[str]
    variants: list              # list[VariantRecord]
    not_in_db: List[str]        # genes with no coordinates anywhere
    fetched_online: List[str]   # genes resolved via mygene.info

@dataclass
class TextEvent:
    text: str

@dataclass
class ErrorEvent:
    message: str

ChatEvent = Union[DatabaseLookupEvent, DatabaseResultEvent,
                  ToolCallEvent, ToolResultEvent, TextEvent, ErrorEvent]

# ── Agentic loop ───────────────────────────────────────────────────────────────

MAX_TOOL_ROUNDS = 6        # safety cap on consecutive tool calls per turn
MAX_VARIANTS_IN_RESULT = 50  # variants sent to Claude per tool call (rest summarised)
# Rough token budget: keep history well under 900k so system prompt + response fit.
# 1 token ≈ 4 chars of JSON.  900_000 * 4 = 3_600_000 chars.
HISTORY_CHAR_BUDGET = 3_000_000

# Fields the Anthropic Messages API accepts for each assistant content block type.
# Anything outside this set (parsed_output, citations, caller, …) causes a 400.
_API_BLOCK_FIELDS: dict[str, set[str]] = {
    "text":              {"type", "text"},
    "tool_use":          {"type", "id", "name", "input"},
    "thinking":          {"type", "thinking", "signature"},
    "redacted_thinking": {"type", "data"},
}


def _clean_block(block) -> dict:
    """
    Convert one SDK content block (Pydantic model or dict) to an API-safe dict.
    Uses an allow-list per block type so SDK-internal fields like
    parsed_output (ParsedTextBlock), citations (TextBlock), and
    caller (ToolUseBlock) are silently dropped.
    """
    if hasattr(block, "model_dump"):
        d = block.model_dump()
    elif isinstance(block, dict):
        d = dict(block)
    else:
        return {"type": "text", "text": str(block)}

    allowed = _API_BLOCK_FIELDS.get(d.get("type", ""))
    if allowed:
        d = {k: v for k, v in d.items() if k in allowed}
    return d


def _strip_thinking_from_old_turns(messages: list[dict]) -> None:
    """
    Remove thinking/redacted_thinking blocks from all assistant turns except
    the last one.  Thinking blocks can be hundreds of tokens each and provide
    no value once a turn is complete — only the current turn benefits from
    extended reasoning context.
    Modifies messages in-place.
    """
    thinking_types = {"thinking", "redacted_thinking"}
    # Find assistant message indices
    assistant_indices = [
        i for i, m in enumerate(messages) if m.get("role") == "assistant"
    ]
    # Keep the last assistant turn intact; strip thinking from all earlier ones
    for i in assistant_indices[:-1]:
        content = messages[i].get("content", [])
        if isinstance(content, list):
            stripped = [b for b in content if b.get("type") not in thinking_types]
            if len(stripped) != len(content):
                messages[i] = {**messages[i], "content": stripped}


def _prune_history_if_needed(messages: list[dict]) -> None:
    """
    If the serialised history exceeds HISTORY_CHAR_BUDGET, drop the oldest
    user+assistant turn pairs until it fits.  Always keeps at least the
    most recent user message so the current turn is never lost.
    Modifies messages in-place.
    """
    import json as _json

    def _size() -> int:
        try:
            return len(_json.dumps(messages, default=str))
        except Exception:
            return 0

    while _size() > HISTORY_CHAR_BUDGET and len(messages) > 2:
        # Drop messages[0] and messages[1] (oldest user+assistant pair).
        # But never drop the very last user message (current turn).
        if len(messages) <= 2:
            break
        # Find the first user/assistant pair to drop
        messages.pop(0)
        # After removing the first message, if it was a user turn the next
        # is probably the assistant turn — remove that too.
        if messages and messages[0].get("role") == "assistant":
            messages.pop(0)


def run_chat_turn(
    client: anthropic.Anthropic,
    parser: VCFParser,
    messages: list[dict],
) -> Generator[ChatEvent, None, None]:
    """
    Execute one user turn, handling any number of tool-use rounds.
    Yields ChatEvent objects for the UI to consume.
    `messages` is modified in-place to maintain conversation history.
    """
    from gene_panels import GENE_COORDINATES

    # Retroactively sanitize any stale Pydantic objects already in the history
    # (messages from before this fix or loaded from an old saved profile).
    for i, msg in enumerate(messages):
        content = msg.get("content", "")
        if isinstance(content, list) and any(hasattr(b, "model_dump") for b in content):
            messages[i] = {**msg, "content": [_clean_block(b) for b in content]}

    for _round in range(MAX_TOOL_ROUNDS):
        # Keep history lean before every API call
        _strip_thinking_from_old_turns(messages)
        _prune_history_if_needed(messages)
        try:
            with client.beta.messages.stream(
                model="claude-opus-4-7",
                max_tokens=6000,
                thinking={"type": "adaptive"},
                system=[
                    {
                        "type": "text",
                        "text": CHAT_SYSTEM_PROMPT,
                        "cache_control": {"type": "ephemeral"},
                    }
                ],
                tools=[QUERY_DATABASES_TOOL, QUERY_VCF_TOOL],
                messages=messages,
                betas=["compact-2026-01-12"],   # auto-summarise history near context limit
            ) as stream:
                response = stream.get_final_message()
        except Exception as e:
            yield ErrorEvent(message=str(e))
            return

        # Append assistant response to history using the allow-list serializer.
        # Compaction blocks (type="server_tool_use" / "server_tool_result") must be
        # preserved verbatim so the server can replace compacted history on the next turn.
        cleaned = []
        for b in response.content:
            d = b.model_dump() if hasattr(b, "model_dump") else dict(b)
            block_type = d.get("type", "")
            allowed = _API_BLOCK_FIELDS.get(block_type)
            if allowed:
                d = {k: v for k, v in d.items() if k in allowed}
            # Blocks not in our allow-list (compaction blocks, etc.) are kept as-is
            cleaned.append(d)
        messages.append({"role": "assistant", "content": cleaned})

        # Separate tool-use blocks from text blocks
        tool_uses = [b for b in response.content if b.type == "tool_use"]
        text_blocks = [b for b in response.content if b.type == "text" and b.text.strip()]

        if not tool_uses:
            # No more tool calls — yield final text and finish
            for block in text_blocks:
                yield TextEvent(text=block.text)
            return

        # ── Execute each tool call ─────────────────────────────────────────
        tool_results = []
        for tu in tool_uses:
            from gene_panels import GENE_COORDINATES as _COORDS

            reason: str = tu.input.get("reason", "")

            # ── Branch: query_databases ────────────────────────────────────
            if tu.name == "query_databases":
                db_genes: list[str] = [g.upper() for g in tu.input.get("genes", [])]
                db_list: list[str] = tu.input.get("databases") or ["clinvar", "gwas", "pharmgkb"]

                yield DatabaseLookupEvent(
                    genes=db_genes, databases=db_list, reason=reason
                )

                try:
                    db_results = lookup_gene_databases(db_genes, databases=db_list)
                except Exception as exc:
                    db_results = {}
                    yield ErrorEvent(message=f"Database lookup error: {exc}")

                # Summarise counts and collect all rsIDs for the UI
                counts: dict = {}
                all_rsids: list[str] = []
                for gene, sources in db_results.items():
                    counts[gene] = {}
                    for src, items in sources.items():
                        if isinstance(items, list):
                            counts[gene][src] = len(items)
                            for item in items:
                                rs = item.get("rsid")
                                if rs and rs not in all_rsids:
                                    all_rsids.append(rs)

                yield DatabaseResultEvent(
                    genes=db_genes, databases=db_list,
                    counts=counts, rsids_found=all_rsids,
                )

                tool_results.append({
                    "type": "tool_result",
                    "tool_use_id": tu.id,
                    "content": json.dumps(
                        {
                            "databases_queried": db_list,
                            "genes_queried": db_genes,
                            "reason": reason,
                            "results": db_results,
                            "total_rsids_found": len(all_rsids),
                            "next_step": (
                                "Use the rsIDs returned above as input to query_vcf "
                                "in rsID mode to look up the patient's genotype at those positions."
                                if all_rsids else
                                "No rsIDs found in these databases for the queried genes. "
                                "Proceed with query_vcf using gene names."
                            ),
                        },
                        default=str,
                    ),
                })
                continue   # move to next tool use block

            # ── Branch: query_vcf ──────────────────────────────────────────
            genes: list[str] = [g.upper() for g in tu.input.get("genes", [])]
            rsids: list[str] = tu.input.get("rsids", []) or []
            include_common: bool = bool(tu.input.get("include_common", False))

            # Label for UI display
            label_parts = []
            if genes:
                label_parts.append(", ".join(genes))
            if rsids:
                label_parts.append(f"rsIDs: {', '.join(rsids)}")
            if include_common:
                label_parts.append("(unfiltered)")
            display_genes = genes + rsids

            # Pre-resolve gene coordinates to detect online lookups
            coord_map = get_coordinates_batch(genes) if genes else {}
            in_hardcoded = set(_COORDS.keys())
            fetched_online = [
                g for g, c in coord_map.items()
                if c is not None and g not in in_hardcoded
            ]
            no_coords = [g for g, c in coord_map.items() if c is None]

            yield ToolCallEvent(
                genes=display_genes, reason=reason, fetched_online=fetched_online
            )

            variants: list[VariantRecord] = []
            try:
                if rsids:
                    # Mode 3: direct rsID lookup (no filter)
                    variants = parser.get_variants_by_rsids(rsids)
                    if genes:
                        # Also do gene query and merge
                        if include_common:
                            gene_vars = parser.get_unfiltered_variants(genes)
                        else:
                            gene_vars = parser.get_variants_for_question(genes)
                        seen_keys = {(v.chrom, v.pos, v.ref, v.alt) for v in variants}
                        for v in gene_vars:
                            if (v.chrom, v.pos, v.ref, v.alt) not in seen_keys:
                                variants.append(v)
                elif include_common:
                    # Mode 2: all variants, no filter
                    variants = parser.get_unfiltered_variants(genes)
                else:
                    # Mode 1: clinical filter (default)
                    variants = parser.get_variants_for_question(genes)
            except Exception as e:
                yield ErrorEvent(message=f"VCF query error: {e}")

            not_in_db = no_coords if not parser.has_annotations else []
            yield ToolResultEvent(
                genes=display_genes,
                variants=variants,
                not_in_db=not_in_db,
                fetched_online=fetched_online,
            )

            mode_note = (
                "rsID direct lookup — no frequency/impact filter applied."
                if rsids else
                "Unfiltered — all variants including common SNPs returned."
                if include_common else
                "Clinical filter applied (ClinVar P/LP, HIGH impact, rare/ultra-rare only)."
            )

            # Cap variants sent to Claude to avoid context overflow.
            # Unfiltered queries can return hundreds of common SNPs; we send
            # the top MAX_VARIANTS_IN_RESULT and summarise the rest.
            variants_for_claude = variants[:MAX_VARIANTS_IN_RESULT]
            truncated = len(variants) - len(variants_for_claude)

            result_payload = {
                "query_mode": "rsid" if rsids else ("unfiltered" if include_common else "clinical"),
                "mode_note": mode_note,
                "genes_queried": genes,
                "rsids_queried": rsids,
                "reason": reason,
                "variants_found": len(variants),
                "variants_shown": len(variants_for_claude),
                "variants_truncated": truncated,
                "genes_resolved_online": fetched_online,
                "genes_without_any_coordinates": no_coords,
                "annotation_scan_available": parser.has_annotations,
                "variants": [v.to_dict() for v in variants_for_claude],
                "interpretation_note": (
                    f"Showing top {len(variants_for_claude)} of {len(variants)} variants "
                    f"({truncated} additional variants omitted to stay within context limits). "
                    "Prioritised by clinical significance."
                    if truncated else
                    "No variants found at these positions/genes — patient likely "
                    "carries the reference genotype."
                    if not variants else ""
                ),
            }

            tool_results.append(
                {
                    "type": "tool_result",
                    "tool_use_id": tu.id,
                    "content": json.dumps(result_payload, default=str),
                }
            )

        messages.append({"role": "user", "content": tool_results})

    # Exceeded MAX_TOOL_ROUNDS
    yield ErrorEvent(message="Maximum tool-call rounds reached. Please rephrase your question.")
