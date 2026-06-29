"""
User profile management.

Stores per-user VCF paths, panel results, chat history, and settings in:
  ~/.genome_vcf_evaluator/
    config.json                  ← global settings (API key)
    gene_coords_cache.json       ← gene lookup cache (shared)
    profiles/
      {name}.json                ← per-user data
      {name}/                    ← per-user VCF storage
        genome.vcf.gz
        genome.vcf.gz.tbi
"""

from __future__ import annotations

import json
import subprocess
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

BASE_DIR = Path.home() / ".genome_vcf_evaluator"
PROFILES_DIR = BASE_DIR / "profiles"
CONFIG_FILE = BASE_DIR / "config.json"


# ── Directory setup ────────────────────────────────────────────────────────────

def _ensure_dirs() -> None:
    BASE_DIR.mkdir(parents=True, exist_ok=True)
    PROFILES_DIR.mkdir(parents=True, exist_ok=True)


def profile_vcf_dir(name: str) -> Path:
    d = PROFILES_DIR / name
    d.mkdir(parents=True, exist_ok=True)
    return d


# ── Profile CRUD ───────────────────────────────────────────────────────────────

def list_profiles() -> list[str]:
    _ensure_dirs()
    return sorted(p.stem for p in PROFILES_DIR.glob("*.json"))


def load_profile(name: str) -> dict:
    _ensure_dirs()
    path = PROFILES_DIR / f"{name}.json"
    if path.exists():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            pass
    return _empty_profile(name)


def save_profile(name: str, data: dict) -> None:
    _ensure_dirs()
    data["name"] = name
    data["updated"] = datetime.now().isoformat()
    path = PROFILES_DIR / f"{name}.json"
    path.write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")


def delete_profile(name: str) -> None:
    _ensure_dirs()
    json_path = PROFILES_DIR / f"{name}.json"
    vcf_dir = PROFILES_DIR / name
    if json_path.exists():
        json_path.unlink()
    if vcf_dir.exists():
        shutil.rmtree(vcf_dir, ignore_errors=True)


def _empty_profile(name: str) -> dict:
    return {
        "name": name,
        "vcf_path": None,
        "created": datetime.now().isoformat(),
        "updated": datetime.now().isoformat(),
        "panel_results": {},     # panel_name → {"variants": [...], "analysis": str, "date": str}
        "chat_api_messages": [], # full Claude message history (serialized)
        "chat_display": [],      # UI-friendly display history
    }


# ── Global config (API key, shared settings) ───────────────────────────────────

def load_config() -> dict:
    _ensure_dirs()
    if CONFIG_FILE.exists():
        try:
            return json.loads(CONFIG_FILE.read_text(encoding="utf-8"))
        except Exception:
            pass
    return {}


def save_config(config: dict) -> None:
    _ensure_dirs()
    CONFIG_FILE.write_text(json.dumps(config, indent=2), encoding="utf-8")


def get_api_key() -> str:
    return load_config().get("api_key", "")


def save_api_key(key: str) -> None:
    config = load_config()
    config["api_key"] = key.strip()
    save_config(config)


# ── VCF file management ────────────────────────────────────────────────────────

def save_vcf_for_profile(
    name: str,
    vcf_bytes: bytes,
    filename: str,
) -> tuple[str, bool]:
    """
    Write uploaded VCF bytes to the profile's storage directory.
    Attempts to create a tabix index automatically.
    Returns (vcf_path, tbi_created).
    """
    dest_dir = profile_vcf_dir(name)
    vcf_path = dest_dir / filename
    vcf_path.write_bytes(vcf_bytes)
    tbi_created = create_tabix_index(str(vcf_path))
    return str(vcf_path), tbi_created


def create_tabix_index(vcf_path: str) -> bool:
    """
    Run `tabix -p vcf <file>` to create a .tbi index.
    Returns True if successful.
    Requires htslib to be installed (brew install htslib / apt install tabix).
    """
    if not vcf_path.endswith(".gz"):
        return False
    tbi_path = vcf_path + ".tbi"
    if Path(tbi_path).exists():
        return True  # Already indexed
    try:
        result = subprocess.run(
            ["tabix", "-p", "vcf", vcf_path],
            capture_output=True,
            timeout=600,   # 10 min — WGS indexing can take a while
        )
        return result.returncode == 0 and Path(tbi_path).exists()
    except FileNotFoundError:
        return False   # tabix not installed
    except subprocess.TimeoutExpired:
        return False


def tabix_available() -> bool:
    """Check if tabix is installed on this system."""
    try:
        subprocess.run(["tabix", "--version"], capture_output=True, timeout=5)
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


# ── Serialization helpers ──────────────────────────────────────────────────────

def serialize_chat_display(display: list[dict]) -> list[dict]:
    """Convert VariantRecord objects → plain dicts for JSON storage."""
    result = []
    for entry in display:
        e = {k: v for k, v in entry.items() if k != "tool_calls"}
        tool_calls = []
        for tc in entry.get("tool_calls", []):
            tc2 = dict(tc)
            tc2["variants"] = [
                v.to_dict() if hasattr(v, "to_dict") else v
                for v in tc2.get("variants", [])
            ]
            tool_calls.append(tc2)
        e["tool_calls"] = tool_calls
        result.append(e)
    return result


def serialize_api_messages(messages: list[dict]) -> list[dict]:
    """Convert Pydantic SDK content blocks → plain dicts for JSON storage."""
    result = []
    for msg in messages:
        m: dict = {"role": msg["role"]}
        content = msg.get("content", "")
        if isinstance(content, str):
            m["content"] = content
        elif isinstance(content, list):
            blocks = []
            for block in content:
                if hasattr(block, "model_dump"):
                    d = block.model_dump()
                elif isinstance(block, dict):
                    d = dict(block)
                else:
                    d = {"type": "text", "text": str(block)}
                # Strip SDK-internal fields (parsed_output, citations, caller…)
                # that the API rejects when the block is echoed back in history.
                block_type = d.get("type", "")
                allowed = {
                    "text":              {"type", "text"},
                    "tool_use":          {"type", "id", "name", "input"},
                    "thinking":          {"type", "thinking", "signature"},
                    "redacted_thinking": {"type", "data"},
                }.get(block_type)
                if allowed:
                    d = {k: v for k, v in d.items() if k in allowed}
                blocks.append(d)
            m["content"] = blocks
        else:
            m["content"] = str(content)
        result.append(m)
    return result


def serialize_panel_variants(variants: list) -> list[dict]:
    """Convert VariantRecord list to JSON-serializable dicts."""
    return [v.to_dict() if hasattr(v, "to_dict") else v for v in variants]
