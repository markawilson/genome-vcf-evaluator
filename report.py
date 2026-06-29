"""
Report generation: markdown output and Rich console display.
"""

from __future__ import annotations

import datetime
from pathlib import Path

from rich.console import Console
from rich.markdown import Markdown
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.rule import Rule
from rich.table import Table
from rich.text import Text

from vcf_parser import VariantRecord

console = Console()


def print_header(vcf_path: str, panels: list[str]) -> None:
    console.print()
    console.print(
        Panel.fit(
            "[bold cyan]Genomic VCF Analyzer[/bold cyan]\n"
            "[dim]Powered by Claude AI — Research & Educational Use Only[/dim]",
            border_style="cyan",
        )
    )
    console.print(f"  [bold]VCF file:[/bold] {vcf_path}")
    console.print(f"  [bold]Panels:[/bold]  {', '.join(panels)}")
    console.print(f"  [bold]Date:[/bold]    {datetime.date.today().isoformat()}")
    console.print()


def print_disclaimer() -> None:
    console.print(
        Panel(
            "[bold yellow]DISCLAIMER[/bold yellow]\n"
            "This analysis is for [bold]research and educational purposes only[/bold]. "
            "It does not constitute medical advice. Clinical decisions should be made by "
            "qualified healthcare professionals using validated clinical testing methods. "
            "Variants of clinical significance should be confirmed by a CLIA-certified "
            "laboratory before any clinical action is taken.",
            border_style="yellow",
        )
    )
    console.print()


def print_variants_table(variants: list[VariantRecord], panel_name: str) -> None:
    if not variants:
        console.print(f"  [dim]No clinically relevant variants found in {panel_name} panel genes.[/dim]")
        return

    table = Table(
        title=f"{panel_name.title()} Panel — {len(variants)} variant(s)",
        show_header=True,
        header_style="bold magenta",
        box=None,
        padding=(0, 1),
    )
    table.add_column("Gene", style="bold", no_wrap=True)
    table.add_column("Variant", style="cyan")
    table.add_column("Zygosity")
    table.add_column("Impact", justify="center")
    table.add_column("ClinVar")
    table.add_column("gnomAD AF", justify="right")

    def impact_style(impact: str) -> str:
        return {
            "HIGH": "bold red",
            "MODERATE": "yellow",
            "LOW": "green",
            "MODIFIER": "dim",
        }.get(impact, "white")

    def clnsig_style(sig: str) -> str:
        s = sig.lower()
        if "pathogenic" in s and "likely" not in s and "conflicting" not in s:
            return "bold red"
        if "likely_pathogenic" in s or "likely pathogenic" in s:
            return "red"
        if "uncertain" in s or "vus" in s:
            return "yellow"
        if "benign" in s:
            return "green"
        return "white"

    for v in variants:
        label = v.hgvs_p or v.hgvs_c or f"{v.ref}>{v.alt}"
        gnomad = f"{v.gnomad_af:.5f}" if v.gnomad_af is not None else "—"
        clnsig = v.clinvar_significance or "—"
        table.add_row(
            v.gene,
            label,
            v.zygosity,
            Text(v.impact, style=impact_style(v.impact)),
            Text(clnsig, style=clnsig_style(clnsig)),
            gnomad,
        )

    console.print(table)
    console.print()


def stream_analysis_to_console(text_iter, panel_name: str) -> str:
    """Stream Claude's analysis to the console, returning the full text."""
    console.print(Rule(f"[bold cyan]{panel_name.title()} Panel Analysis[/bold cyan]"))
    console.print()

    full_text = ""
    # Stream raw markdown — print character by character to simulate streaming
    # Rich can't incrementally render markdown, so we print chunks directly
    for chunk in text_iter:
        console.print(chunk, end="", highlight=False, markup=False)
        full_text += chunk

    console.print()
    console.print()
    return full_text


def save_markdown_report(
    vcf_path: str,
    panels_analyzed: list[str],
    panel_variants: dict[str, list[VariantRecord]],
    panel_analyses: dict[str, str],
    output_path: str,
) -> None:
    lines = [
        "# Genomic VCF Analysis Report",
        "",
        f"**Date**: {datetime.date.today().isoformat()}  ",
        f"**VCF file**: {vcf_path}  ",
        f"**Panels analyzed**: {', '.join(panels_analyzed)}  ",
        "",
        "> **DISCLAIMER**: This analysis is for research and educational purposes only. "
        "It does not constitute medical advice. Clinical decisions should be made by "
        "qualified healthcare professionals using validated clinical testing methods.",
        "",
        "---",
        "",
    ]

    for panel in panels_analyzed:
        variants = panel_variants.get(panel, [])
        analysis = panel_analyses.get(panel, "")

        lines += [
            f"## {panel.title()} Panel",
            "",
            f"**Variants evaluated**: {len(variants)}",
            "",
        ]

        if variants:
            lines += [
                "### Variants Detected",
                "",
                "| Gene | Variant | Zygosity | Impact | ClinVar | gnomAD AF |",
                "|------|---------|----------|--------|---------|-----------|",
            ]
            for v in variants:
                label = v.hgvs_p or v.hgvs_c or f"{v.ref}>{v.alt}"
                gnomad = f"{v.gnomad_af:.6f}" if v.gnomad_af is not None else "—"
                clnsig = v.clinvar_significance or "—"
                lines.append(
                    f"| {v.gene} | {label} | {v.zygosity} | {v.impact} | {clnsig} | {gnomad} |"
                )
            lines.append("")

        if analysis:
            lines += ["### Clinical Interpretation", "", analysis, ""]

        lines += ["---", ""]

    lines += [
        "## About This Report",
        "",
        "This report was generated using the Genomic VCF Analyzer, "
        "which uses Claude (Anthropic) to interpret genomic variants. "
        "Gene panels cover clinically established risk genes. "
        "Variant filtering prioritizes ClinVar-annotated, high-impact, and rare variants.",
        "",
        "*Generated by Genomic VCF Analyzer — Research Use Only*",
    ]

    Path(output_path).write_text("\n".join(lines), encoding="utf-8")
    console.print(f"\n[bold green]Report saved:[/bold green] {output_path}")


def make_spinner(label: str) -> Progress:
    return Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
        transient=True,
    )
