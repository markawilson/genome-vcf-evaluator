#!/usr/bin/env python3
"""
Genomic VCF Analyzer — CLI entry point.

Usage:
  python main.py analyze --vcf sample.vcf.gz --panels cancer cardiovascular
  python main.py query   --vcf sample.vcf.gz --question "What is my APOE genotype?"
"""

import os
import sys

import click
from rich.console import Console

from claude_analyzer import GenomeAnalyzer
from gene_panels import ALL_PANELS
from report import (
    print_header,
    print_disclaimer,
    print_variants_table,
    save_markdown_report,
    stream_analysis_to_console,
    make_spinner,
    console,
)
from vcf_parser import VCFParser

VALID_PANELS = list(ALL_PANELS.keys())


@click.group()
def cli() -> None:
    """Genomic VCF Analyzer — evaluate whole-genome VCF files with Claude AI."""


@cli.command()
@click.option(
    "--vcf",
    required=True,
    type=click.Path(exists=True),
    help="Path to VCF or VCF.gz file.",
)
@click.option(
    "--panels",
    multiple=True,
    type=click.Choice(VALID_PANELS, case_sensitive=False),
    default=VALID_PANELS,
    show_default=True,
    help="Gene panels to analyze. Specify one or more. Defaults to all.",
)
@click.option(
    "--output",
    "-o",
    default=None,
    help="Save markdown report to this path (optional).",
)
@click.option(
    "--api-key",
    envvar="ANTHROPIC_API_KEY",
    default=None,
    help="Anthropic API key (defaults to ANTHROPIC_API_KEY env var).",
)
@click.option(
    "--no-disclaimer",
    is_flag=True,
    default=False,
    help="Suppress the disclaimer banner.",
)
def analyze(
    vcf: str,
    panels: tuple[str, ...],
    output: str | None,
    api_key: str | None,
    no_disclaimer: bool,
) -> None:
    """
    Analyze a VCF file against clinical gene panels.

    Filters variants to clinically relevant findings, then asks Claude to
    interpret them with ACMG/AMP classification and gene/phenotype summaries.
    """
    if not api_key and not os.environ.get("ANTHROPIC_API_KEY"):
        console.print("[bold red]Error:[/bold red] ANTHROPIC_API_KEY is not set.")
        console.print("Set it via --api-key or the ANTHROPIC_API_KEY environment variable.")
        sys.exit(1)

    panels = panels or tuple(VALID_PANELS)
    print_header(vcf, list(panels))
    if not no_disclaimer:
        print_disclaimer()

    try:
        parser = VCFParser(vcf)
        analyzer = GenomeAnalyzer(api_key=api_key)
    except Exception as e:
        console.print(f"[bold red]Initialization error:[/bold red] {e}")
        sys.exit(1)

    if not parser._has_tabix:
        console.print(
            "[yellow]Warning:[/yellow] No tabix index (.tbi) found. "
            "For WGS files, indexing with 'tabix -p vcf sample.vcf.gz' is strongly recommended.\n"
        )

    panel_variants: dict[str, list] = {}
    panel_analyses: dict[str, str] = {}

    for panel_name in panels:
        gene_set = ALL_PANELS[panel_name]
        console.print(f"[bold]Scanning {panel_name} panel[/bold] ({len(gene_set)} genes)…")

        with make_spinner(f"Filtering {panel_name} variants") as progress:
            task = progress.add_task(f"Filtering {panel_name} variants", total=None)
            try:
                variants = parser.filter_panel_variants(panel_name, gene_set)
            except Exception as e:
                console.print(f"[red]Error parsing {panel_name} panel:[/red] {e}")
                panel_variants[panel_name] = []
                panel_analyses[panel_name] = ""
                continue
            progress.update(task, completed=1)

        panel_variants[panel_name] = variants
        print_variants_table(variants, panel_name)

        console.print(f"[dim]Sending {len(variants)} variant(s) to Claude for {panel_name} analysis…[/dim]")
        console.print()

        try:
            full_text = stream_analysis_to_console(
                analyzer.analyze_panel(variants, panel_name, gene_set),
                panel_name,
            )
            panel_analyses[panel_name] = full_text
        except Exception as e:
            console.print(f"[red]Claude API error for {panel_name} panel:[/red] {e}")
            panel_analyses[panel_name] = ""

    if output:
        try:
            save_markdown_report(vcf, list(panels), panel_variants, panel_analyses, output)
        except Exception as e:
            console.print(f"[red]Error saving report:[/red] {e}")

    console.print("\n[bold green]Analysis complete.[/bold green]")


@cli.command()
@click.option(
    "--vcf",
    required=True,
    type=click.Path(exists=True),
    help="Path to VCF or VCF.gz file.",
)
@click.option(
    "--question",
    "-q",
    required=True,
    help='Patient-specific question, e.g. "What is my APOE genotype?"',
)
@click.option(
    "--output",
    "-o",
    default=None,
    help="Save the answer to this file (optional).",
)
@click.option(
    "--api-key",
    envvar="ANTHROPIC_API_KEY",
    default=None,
    help="Anthropic API key (defaults to ANTHROPIC_API_KEY env var).",
)
def query(
    vcf: str,
    question: str,
    output: str | None,
    api_key: str | None,
) -> None:
    """
    Answer a patient-specific question about variants in the VCF.

    Claude will identify the relevant genes, retrieve the variants,
    and provide a targeted clinical answer.

    Examples:
      python main.py query --vcf genome.vcf.gz --question "What is my APOE genotype?"
      python main.py query --vcf genome.vcf.gz --question "Do I carry any BRCA mutations?"
      python main.py query --vcf genome.vcf.gz --question "What is my factor V Leiden status?"
    """
    if not api_key and not os.environ.get("ANTHROPIC_API_KEY"):
        console.print("[bold red]Error:[/bold red] ANTHROPIC_API_KEY is not set.")
        sys.exit(1)

    print_disclaimer()
    console.print(f"[bold]Question:[/bold] {question}\n")

    try:
        analyzer = GenomeAnalyzer(api_key=api_key)
        parser = VCFParser(vcf)
    except Exception as e:
        console.print(f"[bold red]Initialization error:[/bold red] {e}")
        sys.exit(1)

    # Step 1: Ask Claude which genes/rsIDs to look up
    console.print("[dim]Identifying relevant genes for your question…[/dim]")
    try:
        genes_to_query = analyzer.identify_relevant_genes(question)
    except Exception as e:
        console.print(f"[red]Error identifying genes:[/red] {e}")
        genes_to_query = []

    if not genes_to_query:
        console.print("[yellow]Could not identify specific genes. Please mention gene names explicitly.[/yellow]")
        sys.exit(1)

    from gene_panels import GENE_COORDINATES
    coord_genes = [g for g in genes_to_query if g in GENE_COORDINATES]
    annot_genes = [g for g in genes_to_query if g not in GENE_COORDINATES]

    if coord_genes:
        console.print(f"[dim]Fast lookup (coordinates): {', '.join(coord_genes)}[/dim]")
    if annot_genes:
        if parser.has_annotations:
            console.print(f"[dim]Annotation scan: {', '.join(annot_genes)} (may be slower for large files)[/dim]")
        else:
            console.print(
                f"[yellow]Note:[/yellow] Gene(s) {', '.join(annot_genes)} are not in the built-in panel "
                "database and this VCF has no VEP/SnpEff annotations — skipping those genes."
            )
            annot_genes = []

    all_genes = coord_genes + annot_genes
    if not all_genes:
        console.print("[red]No queryable genes found.[/red] Try asking about a specific gene by name.")
        sys.exit(1)

    console.print()

    # Step 2: Retrieve variants from VCF for those genes
    with make_spinner("Retrieving variants from VCF") as progress:
        task = progress.add_task("Retrieving variants", total=None)
        try:
            variants = parser.get_variants_for_question(all_genes)
        except Exception as e:
            console.print(f"[red]Error reading VCF:[/red] {e}")
            sys.exit(1)
        progress.update(task, completed=1)

    if variants:
        print_variants_table(variants, "Query")
    else:
        console.print("[dim]No variants detected in these genes (may indicate reference genotype).[/dim]\n")

    # Step 3: Stream Claude's answer
    console.print("[bold cyan]── Claude's Answer ──[/bold cyan]\n")
    full_answer = ""
    try:
        for chunk in analyzer.answer_question(variants, question, context_genes=all_genes):
            console.print(chunk, end="", highlight=False, markup=False)
            full_answer += chunk
    except Exception as e:
        console.print(f"\n[red]Claude API error:[/red] {e}")
        sys.exit(1)

    console.print("\n")

    if output and full_answer:
        from pathlib import Path
        Path(output).write_text(
            f"# Genomic Query Answer\n\n**Question**: {question}\n\n{full_answer}\n\n"
            f"*Generated by Genomic VCF Analyzer — Research Use Only*\n",
            encoding="utf-8",
        )
        console.print(f"[bold green]Answer saved:[/bold green] {output}")


if __name__ == "__main__":
    cli()
