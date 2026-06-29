"""
Claude API integration for genome variant analysis.
Uses claude-opus-4-7 with adaptive thinking, streaming, and prompt caching.
"""

from __future__ import annotations

import json
import os
from typing import Iterator

import anthropic

from vcf_parser import VariantRecord

SYSTEM_PROMPT = """\
You are an expert clinical geneticist and genomic medicine specialist with deep knowledge of:
- ACMG/AMP variant classification guidelines (Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign)
- ClinVar, gnomAD, OMIM, and other clinical genomic databases
- Gene function, disease mechanisms, and phenotype-genotype correlations
- Pharmacogenomics and personalized medicine
- Population genetics and variant frequency interpretation

When analyzing genomic variants, you will:
1. Classify variants according to ACMG/AMP criteria when possible
2. Explain the clinical significance of each variant clearly
3. Describe the gene function and its role in the relevant disease pathway
4. Note zygosity (heterozygous vs homozygous) and its implications
5. Highlight the most clinically actionable findings prominently
6. Discuss penetrance and expressivity where relevant
7. Mention recommended clinical follow-up when appropriate

Output format:
- Use clear headings for each gene/variant
- Lead with the most clinically significant findings
- Use plain language for clinical significance, but include technical details
- Always note when a variant is of uncertain significance (VUS)

IMPORTANT DISCLAIMER: This analysis is for research and educational purposes only. \
It does not constitute medical advice. Clinical decisions should be made by qualified \
healthcare professionals using validated clinical testing methods. Variants of clinical \
significance should be confirmed by a CLIA-certified laboratory before any clinical action is taken.\
"""

PANEL_DESCRIPTIONS = {
    "cancer": "hereditary cancer risk (BRCA1/2, Lynch syndrome, Li-Fraumeni, and other cancer predisposition genes)",
    "cardiovascular": "cardiovascular genetics (familial hypercholesterolemia, cardiomyopathies, channelopathies, aortopathies)",
    "longevity": "longevity and healthy aging genomics (APOE, FOXO3, telomere biology, mTOR pathway, sirtuins)",
    "mitochondrial": "mitochondrial function and energy metabolism (POLG, respiratory chain, mtDNA-encoded genes)",
    "depression_anxiety": (
        "depression, anxiety, and antidepressant pharmacogenomics — serotonin system (SLC6A4/SERT, HTR2A, HTR1A), "
        "dopamine/catecholamine system (COMT, MAOA, DRD2), norepinephrine system (SLC6A2), "
        "neuroplasticity (BDNF, NTRK2), HPA stress axis (FKBP5, NR3C1, CRHR1), "
        "glutamate/ketamine target (GRIN2B), drug transport (ABCB1/P-gp), "
        "and antidepressant drug metabolism (CYP2D6, CYP2C19)"
    ),
    "male_hormones": (
        "male hormone pathway genomics — androgen signalling (AR, SRD5A2, SRD5A1), "
        "steroidogenesis (STAR, CYP11A1, HSD3B1, CYP17A1, HSD17B3, NR5A1), "
        "oestrogen conversion and binding (CYP19A1/aromatase, ESR1, ESR2), "
        "SHBG and free testosterone, HPG axis (GNRH1, KISS1R, LHCGR, FSHR), "
        "testicular markers (INSL3), and metabolic axis (IGF1, CYP11B1)"
    ),
    "female_hormones": (
        "female hormone pathway genomics — oestrogen signalling (ESR1, ESR2, CYP19A1, HSD17B1, CYP1B1, COMT), "
        "progesterone signalling (PGR), HPG axis (GNRH1, FSHB, FSHR, LHCGR), "
        "ovarian reserve and AMH axis (AMH, AMHR2, GDF9), "
        "adrenal steroidogenesis and androgen excess (CYP17A1, CYP11A1, CYP21A2), "
        "stress-hormone interaction (HSD11B1, NR3C1), prolactin (PRLR), thyroid axis (TSHR), "
        "and thrombosis risk for OCP/HRT counselling (F5 Factor V Leiden, F2 prothrombin)"
    ),
}


class GenomeAnalyzer:
    def __init__(self, api_key: str | None = None) -> None:
        self.client = anthropic.Anthropic(
            api_key=api_key or os.environ.get("ANTHROPIC_API_KEY")
        )

    def _format_variants_for_prompt(self, variants: list[VariantRecord]) -> str:
        if not variants:
            return "No clinically relevant variants detected in this panel's genes."

        lines = [f"Total variants to analyze: {len(variants)}\n"]
        for i, v in enumerate(variants, 1):
            d = v.to_dict()
            fallback = f"{d['ref']}>{d['alt']}"
            label = d["variant"] or fallback
            lines.append(f"### Variant {i}: {d['gene']} — {label}")
            lines.append(f"- **Genomic position**: {d['chrom']}:{d['pos']}")
            if d["rsid"] and d["rsid"] != ".":
                lines.append(f"- **rsID**: {d['rsid']}")
            lines.append(f"- **Zygosity**: {d['zygosity']}")
            if d["clinvar_significance"]:
                lines.append(f"- **ClinVar**: {d['clinvar_significance']}")
            lines.append(f"- **Consequence**: {d['consequence']} ({d['impact']} impact)")
            if d["hgvs_c"]:
                lines.append(f"- **HGVS coding**: {d['hgvs_c']}")
            if d["hgvs_p"]:
                lines.append(f"- **Protein change**: {d['hgvs_p']}")
            if d["gnomad_af"] is not None:
                lines.append(f"- **gnomAD population frequency**: {d['gnomad_af']:.6f}")
            if d["depth"] is not None:
                lines.append(f"- **Sequencing depth**: {d['depth']}x")
            if d["quality"] is not None:
                lines.append(f"- **Variant quality**: {d['quality']:.1f}")
            lines.append(f"- **Filter status**: {d['filter_status']}")
            lines.append("")

        return "\n".join(lines)

    def analyze_panel(
        self,
        variants: list[VariantRecord],
        panel_name: str,
        panel_genes: set[str],
    ) -> Iterator[str]:
        """
        Stream a clinical interpretation of the panel variants.
        Yields text chunks as they arrive from the API.
        """
        panel_desc = PANEL_DESCRIPTIONS.get(panel_name, panel_name)
        panel_display = panel_name.replace("_", " ").title()
        genes_checked = ", ".join(sorted(panel_genes))
        variant_text = self._format_variants_for_prompt(variants)

        # Build panel-specific analysis instructions
        if panel_name == "depression_anxiety":
            extra_instructions = """\

For this panel, also include:
7. **Antidepressant Pharmacogenomics** — For CYP2D6 and CYP2C19 variants found, predict metaboliser \
status (poor/intermediate/normal/ultra-rapid) and list which antidepressants are affected \
(e.g. fluoxetine, paroxetine, escitalopram, sertraline, venlafaxine, amitriptyline). \
For ABCB1 variants, comment on CNS drug penetration. For SLC6A4 and HTR2A variants, \
comment on SSRI response likelihood.
"""
        elif panel_name in ("male_hormones", "female_hormones"):
            extra_instructions = """\

For this panel, also include:
7. **Hormonal Pathway Impact** — Describe how any variants found might affect hormone levels, \
signalling, or metabolism. Comment on free vs. total hormone ratios (SHBG variants), \
aromatase activity (CYP19A1), and receptor sensitivity. For female panel: note any variants \
relevant to contraceptive or HRT safety (F5/F2 thrombosis risk, ESR1 receptor sensitivity).
"""
        elif panel_name == "neurological":
            extra_instructions = """\

For this panel, also include:
7. **APOE Diplotype & Alzheimer's Risk** — Look for rs429358 and rs7412 in the variant list above \
and call the APOE diplotype using this table (columns = rs429358 genotype; rows = rs7412 genotype):

| rs7412 \\ rs429358 | T/T (hom ref) | C/T (het) | C/C (hom alt) |
|--------------------|--------------|-----------|--------------|
| C/C (hom ref)      | ε3/ε3        | ε3/ε4     | ε4/ε4        |
| C/T (het)          | ε2/ε3        | ε2/ε4 *   | ε3/ε4 *      |
| T/T (hom alt)      | ε2/ε2        | ε2/ε3 *   | ε2/ε4 *      |

(* when both sites are heterozygous, the Arg112/Cys158 combination — "ε1" — does not occur naturally, \
so the most parsimonious phased call is given.)

If neither rs429358 nor rs7412 is in the variant list, state "rs429358 and rs7412 not called in this VCF \
— APOE diplotype cannot be determined."

After calling the diplotype, interpret it in the neurological context:
- State the diplotype clearly (e.g. "APOE ε3/ε4")
- Give the approximate lifetime Alzheimer's risk relative to ε3/ε3 baseline \
  (ε2/ε3 ≈ 25% lower; ε3/ε3 = baseline; ε3/ε4 ≈ 2–3× higher; ε4/ε4 ≈ 8–12× higher)
- State the approximate mean age of onset shift (ε4/ε4 approximately one decade earlier)
- Mention the relevance to Lewy body dementia and TBI recovery
- Note the implications for first-degree relatives (autosomal co-dominant)
- Clarify that APOE ε4 is a risk modifier, not a deterministic mutation; many ε4 carriers never develop AD
"""
        elif panel_name == "cardiovascular":
            extra_instructions = """\

For this panel, also include:
7. **APOE Diplotype & Cardiovascular Risk** — Look for rs429358 and rs7412 in the variant list above \
and call the APOE diplotype using this table (columns = rs429358 genotype; rows = rs7412 genotype):

| rs7412 \\ rs429358 | T/T (hom ref) | C/T (het) | C/C (hom alt) |
|--------------------|--------------|-----------|--------------|
| C/C (hom ref)      | ε3/ε3        | ε3/ε4     | ε4/ε4        |
| C/T (het)          | ε2/ε3        | ε2/ε4 *   | ε3/ε4 *      |
| T/T (hom alt)      | ε2/ε2        | ε2/ε3 *   | ε2/ε4 *      |

(* when both sites are heterozygous, the Arg112/Cys158 combination — "ε1" — does not occur naturally, \
so the most parsimonious phased call is given.)

If neither rs429358 nor rs7412 is in the variant list, state "rs429358 and rs7412 not called in this VCF \
— APOE diplotype cannot be determined."

After calling the diplotype, interpret it in the cardiovascular context:
- State the diplotype clearly (e.g. "APOE ε3/ε4")
- Describe the LDL-C effect: ε4 raises LDL ~10–15 mg/dL per allele via impaired LDL-receptor recycling; \
  ε2 lowers LDL but ε2/ε2 homozygotes risk type III hyperlipoproteinemia (dysbetalipoproteinemia — \
  markedly elevated IDL and triglycerides, xanthomas, premature CAD)
- Note the implication for statin therapy: ε4 carriers may have attenuated LDL response to statins \
  and may benefit from more intensive lipid-lowering
- For ε3/ε4 or ε4/ε4: recommend fasting lipid panel with IDL/VLDL fractionation if not already done
- Note that APOE status should be interpreted alongside any LDLR, APOB, or PCSK9 findings in this panel
"""
        else:
            extra_instructions = ""

        user_message = f"""\
## Genomic Panel Analysis: {panel_display}

**Panel scope**: {panel_desc}
**Genes analyzed**: {genes_checked}
**Variants found** (pre-filtered for potential clinical relevance):

{variant_text}

Please provide a comprehensive clinical interpretation of these findings. Structure your response as:

1. **Executive Summary** — Most important findings in 2-3 sentences
2. **Clinically Significant Variants** — Detailed analysis of Pathogenic/Likely Pathogenic variants
3. **Variants of Uncertain Significance (VUS)** — Brief discussion of VUS findings
4. **Gene & Pathway Summary** — Overview of the genes and pathways examined
5. **Recommendations** — Suggested clinical follow-up (genetic counseling, confirmatory testing, surveillance)
6. **Negative Findings** — Notable genes where no concerning variants were detected\
{extra_instructions}
If no variants were found, provide a reassuring summary and explain what the panel covers.\
"""

        with self.client.messages.stream(
            model="claude-opus-4-7",
            max_tokens=8000,
            thinking={"type": "adaptive"},
            system=[
                {
                    "type": "text",
                    "text": SYSTEM_PROMPT,
                    "cache_control": {"type": "ephemeral"},
                }
            ],
            messages=[{"role": "user", "content": user_message}],
        ) as stream:
            for text in stream.text_stream:
                yield text

    def answer_question(
        self,
        variants: list[VariantRecord],
        question: str,
        context_genes: list[str] | None = None,
    ) -> Iterator[str]:
        """
        Stream an answer to a patient-specific genomic question.
        Yields text chunks as they arrive from the API.
        """
        variant_text = self._format_variants_for_prompt(variants)
        gene_context = ""
        if context_genes:
            gene_context = f"\n**Genes queried**: {', '.join(context_genes)}\n"

        user_message = f"""\
## Patient-Specific Genomic Query

**Question**: {question}
{gene_context}
**Relevant variants found in the genome**:

{variant_text}

Please answer the patient's question directly and thoroughly, including:
- The specific genotype/variant information relevant to the question
- Clinical significance and what this means for the individual
- Any important caveats or limitations of this information
- Whether confirmatory clinical testing is recommended

Be specific about the actual variants found (or their absence).\
"""

        with self.client.messages.stream(
            model="claude-opus-4-7",
            max_tokens=4000,
            thinking={"type": "adaptive"},
            system=[
                {
                    "type": "text",
                    "text": SYSTEM_PROMPT,
                    "cache_control": {"type": "ephemeral"},
                }
            ],
            messages=[{"role": "user", "content": user_message}],
        ) as stream:
            for text in stream.text_stream:
                yield text

    def generate_summary(
        self,
        patient_name: str,
        panel_results: dict,
        chat_display: list[dict],
    ):
        """
        Stream a comprehensive patient summary of all findings.
        panel_results: {panel_name: {"analysis": str, "variants": list}}
        chat_display:  the UI display message list from session state
        """
        from datetime import date

        panel_sections = []
        for panel_name, result in panel_results.items():
            analysis = result.get("analysis", "").strip()
            variants = result.get("variants", [])
            if not analysis and not variants:
                continue
            var_summary = (
                f"{len(variants)} clinically filtered variant(s) found"
                if variants else "No clinically relevant variants detected"
            )
            panel_sections.append(
                f"### {panel_name.title()} Panel\n"
                f"*{var_summary}*\n\n"
                f"{analysis[:3000]}"
            )

        chat_qa_pairs = []
        i = 0
        display = list(chat_display)
        while i < len(display):
            entry = display[i]
            if entry.get("role") == "user" and entry.get("text"):
                question_text = entry["text"]
                answer_text = ""
                if i + 1 < len(display) and display[i + 1].get("role") == "assistant":
                    answer_text = display[i + 1].get("text", "")[:800]
                if question_text and answer_text:
                    chat_qa_pairs.append(
                        f"**Q: {question_text}**\n{answer_text}"
                    )
                i += 2
            else:
                i += 1

        panels_text = "\n\n".join(panel_sections) if panel_sections else "No panel analyses completed."
        chat_text = "\n\n---\n\n".join(chat_qa_pairs[:6]) if chat_qa_pairs else "No chat queries completed."

        user_message = f"""\
Generate a comprehensive clinical genomics summary for patient: **{patient_name}**
Date: {date.today().isoformat()}

## Panel Analysis Results

{panels_text}

## Genomic Chat Q&A Highlights

{chat_text}

---

Please produce a well-structured summary with these sections:

1. **Executive Summary** — 3–5 sentences covering the most important findings overall
2. **Key Findings by Domain** — organized by clinical area (cancer risk, cardiovascular, pharmacogenomics, etc.)
3. **Clinically Actionable Variants** — variants warranting follow-up, with recommended action
4. **Reassuring Findings** — important negative results worth noting
5. **Recommendations** — specific suggested next steps (clinical genetics referral, confirmatory testing, surveillance)
6. **Limitations of This Analysis** — what this analysis does and does not cover

Be concise, clear, and clinically useful. Avoid repeating the same information multiple times.\
"""

        with self.client.messages.stream(
            model="claude-opus-4-7",
            max_tokens=5000,
            thinking={"type": "adaptive"},
            system=[
                {
                    "type": "text",
                    "text": SYSTEM_PROMPT,
                    "cache_control": {"type": "ephemeral"},
                }
            ],
            messages=[{"role": "user", "content": user_message}],
        ) as stream:
            for text in stream.text_stream:
                yield text

    def analyze_ancestry(
        self,
        rsid_variants: list[VariantRecord],
        gene_variants: list[VariantRecord],
    ) -> Iterator[str]:
        """
        Stream an ancestry interpretation based on ancestry-informative markers.
        rsid_variants: results from direct rsID lookup (key AIMs).
        gene_variants: unfiltered variants from ancestry gene regions
                       (includes full chrM scan via MT_GENOME entry).
        """
        _MITO_CHROMS = {"chrM", "MT", "M", "chrMT"}

        def _is_mito(v: VariantRecord) -> bool:
            return v.gene == "MT_GENOME" or v.chrom in _MITO_CHROMS

        def _fmt(variants: list[VariantRecord]) -> str:
            if not variants:
                return "None detected."
            lines = []
            for v in variants:
                d = v.to_dict()
                rsid = d["rsid"] if d["rsid"] and d["rsid"] != "." else "—"
                af_str = f"gnomAD AF {d['gnomad_af']:.4f}" if d["gnomad_af"] is not None else ""
                hgvs = d.get("hgvs_p") or d.get("hgvs_c") or f"{d['ref']}>{d['alt']}"
                lines.append(
                    f"- {d['gene']} {hgvs} ({rsid}) | {d['zygosity']} | {af_str}"
                )
            return "\n".join(lines)

        def _fmt_mito(variants: list[VariantRecord]) -> str:
            """
            Format chrM variants with explicit position notation so Claude can
            match them to the mtDNA haplogroup tree (positions relative to rCRS).
            """
            mito = [v for v in variants if _is_mito(v)]
            if not mito:
                return (
                    "No mitochondrial variants detected. The VCF may not include "
                    "chrM/MT calls, or the sample may be haplogroup H (identical to "
                    "the rCRS reference, so no variants would appear)."
                )
            # Sort by position so the pattern reads naturally along the chromosome
            mito.sort(key=lambda v: v.pos)
            lines = []
            for v in mito:
                d = v.to_dict()
                rsid = d["rsid"] if d["rsid"] and d["rsid"] != "." else "—"
                # Prefer annotated HGVS_c (often "m.XXXXREF>ALT" from VEP/SnpEff),
                # fall back to building it from position
                hgvs = d.get("hgvs_c") or f"m.{d['pos']}{d['ref']}>{d['alt']}"
                af_str = f"(gnomAD AF {d['gnomad_af']:.4f})" if d["gnomad_af"] is not None else ""
                lines.append(f"- chrM:{d['pos']:>5}  {hgvs:<28}  ({rsid}) {af_str}")
            return "\n".join(lines)

        # Partition variants: chrM goes to haplogroup section, everything else to AIM section
        all_variants = rsid_variants + gene_variants
        non_mito_gene = [v for v in gene_variants if not _is_mito(v)]
        # rsid_variants may also contain chrM anchor rsIDs — include them in mito text
        mito_text = _fmt_mito(all_variants)

        rsid_text = _fmt(rsid_variants)
        gene_text = _fmt(non_mito_gene[:80])   # cap to avoid context overflow

        user_message = f"""\
## Genomic Ancestry Analysis

Please provide a thorough genomic ancestry interpretation structured as follows.

---

### 1. Maternal Lineage — mtDNA Haplogroup

**Mitochondrial variants found** (full chrM region scan, positions relative to rCRS):
{mito_text}

Assign the most specific maternal haplogroup supportable by the evidence above.

**Haplogroup assignment guide** (rCRS is haplogroup H2a2a1 — haplogroup H individuals \
show few or no variants vs reference; other haplogroups accumulate branch-defining variants):

Key macro-branch positions:
- **Clade M** (Asian/Oceanian): m.10400T, m.14783C, m.15043A
  - Sub-clades: D (m.5178C>A, East Asian), C (m.13263G), G (m.4833G)
- **Clade N** (non-M, non-African): m.10873C
  - **Clade R** (within N): m.12705T
    - **HV** (within R): m.14766C → ancestor of H and V
      - **H** (~45% European): rCRS-like; H1 adds m.3010A; H3 adds m.6776C
      - **V** (~5% European/Saami): m.4580A, m.15904C
    - **J** (~12% European/Middle Eastern): m.295T, m.3010A, m.10398G, m.13708A
    - **T** (~8% European/Middle Eastern): m.4216C, m.11251G, m.13368A
    - **K** (~6% European/Middle Eastern): m.9055A (subclade of U8)
    - **U** (~7% European, pre-Neolithic): m.11467G, m.12308G, m.12372A
      - U5 (oldest European clade): m.150T, m.7768G
    - **W** (~3% European/Central Asian): m.11947G, m.15884C
    - **X** (~2% European/Native American): m.6221T, m.14470C
    - **I** (~2% European): m.8494C, m.15607G
  - **A** (East Asian/Native American): m.663G, m.4248C
  - **B** (East Asian/Native American): m.8281-8289 deletion
  - **F** (East Asian): m.3970T
- **African L clades** (basal): L0 (m.10810T), L1 (m.3516A), L2 (m.16278T), \
  L3 (ancestral to M+N) — many variants vs rCRS

State: (a) the haplogroup call, (b) the population most associated with it and \
approximate frequency, (c) confidence given how many chrM variants are available, \
(d) if chrM is absent or coverage is too sparse to call, say so explicitly.

---

### 2. Continental Ancestry Signals (autosomal AIMs)

**Directly queried ancestry-informative SNPs (rsID lookup)**:
{rsid_text}

**Additional variants from ancestry-informative gene regions (unfiltered)**:
{gene_text}

Discuss what the detected autosomal marker genotypes suggest about continental ancestry \
proportions (European, sub-Saharan African, East Asian, South Asian, \
Native American/Amerindian, Middle Eastern, Oceanian). \
Focus on markers with strong population differentiation (high Fst).

---

### 3. Key Marker Breakdown
For each important autosomal marker found, explain:
- What population(s) the allele tags
- The patient's genotype and zygosity
- Population frequency context (e.g. "this allele is present in ~99% of Europeans \
but <5% of sub-Saharan Africans")

### 4. Notable Absent Markers
Mention key markers that were NOT detected (homozygous reference), and what that tells us.

### 5. Phenotypic Implications
Based on the pigmentation, morphology, and metabolism markers, discuss likely \
phenotypic traits (skin tone, eye/hair colour tendency, lactase status, \
alcohol metabolism).

### 6. Limitations
Be explicit that:
- This is a *genomic ancestry estimate* — not a cultural or ethnic identity statement
- A handful of AIMs gives directional signals, not a precise percentage breakdown
- For ancestry percentages, dedicated tools (e.g. admixture analysis with 1000 Genomes \
reference panel) are required
- Genetic ancestry and self-identified ancestry/ethnicity are distinct concepts
- Absence of a variant at a queried position may mean homozygous reference *or* \
that the position was not called in this VCF
- The mtDNA haplogroup reflects maternal lineage only — one thread of ancestry among many

IMPORTANT DISCLAIMER: This analysis is for research and educational purposes only. \
It does not constitute medical advice or a definitive ancestry determination.\
"""

        with self.client.messages.stream(
            model="claude-opus-4-7",
            max_tokens=5000,
            thinking={"type": "adaptive"},
            system=[
                {
                    "type": "text",
                    "text": SYSTEM_PROMPT,
                    "cache_control": {"type": "ephemeral"},
                }
            ],
            messages=[{"role": "user", "content": user_message}],
        ) as stream:
            for text in stream.text_stream:
                yield text

    def identify_relevant_genes(self, question: str) -> list[str]:
        """
        Use Claude to extract gene names mentioned in a natural-language question.
        Returns a list of gene symbols to look up.
        """
        response = self.client.messages.create(
            model="claude-opus-4-7",
            max_tokens=500,
            messages=[
                {
                    "role": "user",
                    "content": (
                        f"Extract all human gene symbols mentioned in this genomic question. "
                        f"Also include genes commonly associated with the phenotypes/conditions mentioned. "
                        f"Return ONLY a JSON array of gene symbols (uppercase), nothing else.\n\n"
                        f"Question: {question}"
                    ),
                }
            ],
        )
        text = response.content[0].text.strip()
        # Extract JSON array even if surrounded by markdown
        start = text.find("[")
        end = text.rfind("]") + 1
        if start >= 0 and end > start:
            try:
                genes = json.loads(text[start:end])
                return [g.upper() for g in genes if isinstance(g, str)]
            except json.JSONDecodeError:
                pass
        # Fallback: split on commas/spaces
        return [w.strip('[]",\' ').upper() for w in text.split() if w.strip('[]",\' ')]
