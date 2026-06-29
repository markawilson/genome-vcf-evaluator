"""
Genomic VCF Analyzer — Streamlit browser UI.

Run with:
    streamlit run app.py
"""

from __future__ import annotations

import os
from datetime import date
from pathlib import Path
from typing import Optional

import streamlit as st

from chat_engine import (
    ChatEvent, DatabaseLookupEvent, DatabaseResultEvent,
    ToolCallEvent, ToolResultEvent, TextEvent, ErrorEvent,
    run_chat_turn,
)
from claude_analyzer import GenomeAnalyzer
from gene_panels import ALL_PANELS, GENE_COORDINATES, ANCESTRY_GENES, ANCESTRY_RSIDS
from auth import register_user, verify_login, verify_admin, list_users, delete_user, reset_user_password
from profile_manager import (
    list_profiles, load_profile, save_profile, delete_profile,
    save_vcf_for_profile, create_tabix_index, tabix_available,
    get_api_key, save_api_key,
    serialize_chat_display, serialize_api_messages, serialize_panel_variants,
    PROFILES_DIR, user_profiles_dir,
)
from report import save_markdown_report
from vcf_parser import VCFParser, VariantRecord

# ── Page config ────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Genomic VCF Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
.disclaimer-box {
    background:#fff3cd; border:1px solid #ffc107; border-radius:8px;
    padding:12px 16px; margin:8px 0 16px 0; font-size:0.85rem; color:#856404;
}
.profile-badge {
    background:#e8f4f8; border:1px solid #bee3f8; border-radius:6px;
    padding:4px 10px; font-size:0.85rem; color:#2b6cb0; display:inline-block;
}
</style>
""", unsafe_allow_html=True)

# ── Authentication ─────────────────────────────────────────────────────────────

if "logged_in_user" not in st.session_state:
    st.session_state["logged_in_user"] = None

def _show_login_page():
    # Inject login page styles
    st.markdown("""
    <style>
    .login-header {
        text-align: center;
        padding: 2rem 0 0.5rem 0;
    }
    .login-header h1 {
        font-size: 2.2rem;
        margin-bottom: 0.2rem;
    }
    .login-subtitle {
        text-align: center;
        color: #888;
        font-size: 0.95rem;
        margin-bottom: 1.5rem;
    }
    .dna-art {
        text-align: center;
        font-size: 3.5rem;
        letter-spacing: 0.3rem;
        line-height: 1.2;
        padding: 1rem 0;
        opacity: 0.85;
    }
    .helix-row {
        font-family: monospace;
        font-size: 0.85rem;
        color: #4a9;
        text-align: center;
        line-height: 1.4;
        letter-spacing: 0.15rem;
    }
    </style>
    """, unsafe_allow_html=True)

    # Spacer to push content down from top
    st.markdown("<br>", unsafe_allow_html=True)

    # Centered layout using columns
    col_left, col_center, col_right = st.columns([1, 2, 1])

    with col_center:
        # DNA helix art
        st.markdown("""
        <div class="dna-art">🧬</div>
        <div class="helix-row">
        ╭─A══T─╮&nbsp;&nbsp;╭─G══C─╮&nbsp;&nbsp;╭─T══A─╮&nbsp;&nbsp;╭─C══G─╮<br>
        ╰─T══A─╯&nbsp;&nbsp;╰─C══G─╯&nbsp;&nbsp;╰─A══T─╯&nbsp;&nbsp;╰─G══C─╯
        </div>
        """, unsafe_allow_html=True)

        st.markdown('<div class="login-header"><h1>Genome VCF Analyzer</h1></div>', unsafe_allow_html=True)
        st.markdown('<div class="login-subtitle">Personal whole-genome variant analysis powered by AI</div>', unsafe_allow_html=True)

        tab_login, tab_register = st.tabs(["Log in", "Create account"])

        with tab_login:
            with st.form("login_form"):
                username = st.text_input("Username", key="login_user")
                password = st.text_input("Password", type="password", key="login_pass")
                submitted = st.form_submit_button("Log in", use_container_width=True, type="primary")
                if submitted:
                    if verify_login(username, password):
                        st.session_state["logged_in_user"] = username.strip().lower()
                        st.rerun()
                    else:
                        st.error("Invalid username or password.")

        with tab_register:
            with st.form("register_form"):
                new_user = st.text_input("Choose a username", key="reg_user")
                new_pass = st.text_input("Choose a password (min 6 chars)", type="password", key="reg_pass")
                new_pass2 = st.text_input("Confirm password", type="password", key="reg_pass2")
                reg_submitted = st.form_submit_button("Create account", use_container_width=True, type="primary")
                if reg_submitted:
                    if new_pass != new_pass2:
                        st.error("Passwords do not match.")
                    else:
                        ok, msg = register_user(new_user, new_pass)
                        if ok:
                            st.success(msg + " You can now log in.")
                        else:
                            st.error(msg)

        st.markdown("<br>", unsafe_allow_html=True)
        st.caption("For research and educational purposes only. Not medical advice.")


if st.session_state["logged_in_user"] is None:
    _show_login_page()
    st.stop()

# The authenticated user — passed to all profile functions
_AUTH_USER = st.session_state["logged_in_user"]


# ── Constants ──────────────────────────────────────────────────────────────────
PANEL_ICONS = {
    "cancer": "🔬", "cardiovascular": "❤️", "longevity": "⏳",
    "mitochondrial": "⚡", "connective_tissue": "🦴",
    "pharmacogenomics": "💊", "neurological": "🧠",
    "depression_anxiety": "🧘", "male_hormones": "♂️", "female_hormones": "♀️",
}
PANEL_DESCRIPTIONS = {
    "cancer":             "BRCA1/2, TP53, MLH1, MSH2, APC, PTEN + 16 more",
    "cardiovascular":     "LDLR, MYH7, SCN5A, KCNQ1, FBN1, TTN + 19 more",
    "longevity":          "APOE, FOXO3, TERT, SIRT1/3/6, MTOR + 6 more",
    "mitochondrial":      "POLG, TFAM, SURF1 (nuclear) + MT-ND/CYB/CO (mtDNA)",
    "connective_tissue":  "EDS: COL5A1/2, COL1A1/2, TNXB, PLOD1… + Loeys-Dietz",
    "pharmacogenomics":   "CYP2D6/2C19/2C9, TPMT, DPYD, VKORC1, SLCO1B1, G6PD…",
    "neurological":       "Parkinson, Alzheimer, ALS, Huntington, ataxias, epilepsy",
    "depression_anxiety": "SLC6A4, HTR2A, COMT, MAOA, BDNF, FKBP5, CYP2D6/2C19 + more",
    "male_hormones":      "AR, SRD5A2, CYP19A1, ESR1, SHBG, LHCGR, FSHR + more",
    "female_hormones":    "ESR1/2, PGR, CYP19A1, AMH, CYP21A2, FSHR, F5, F2 + more",
}

# ── Session state defaults ─────────────────────────────────────────────────────
def _init_state() -> None:
    defaults: dict = {
        "current_profile": None,
        "vcf_path": None,
        "panel_variants": {},
        "panel_analyses": {},
        "analysis_done": False,
        "query_variants": [],
        "query_answer": "",
        "chat_api_messages": [],
        "chat_display": [],
        "summary_text": "",
        "ancestry_analysis": "",
        "ancestry_rsid_variants": [],
        "ancestry_gene_variants": [],
        "editing_api_key": False,
        "show_vcf_upload": False,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v

_init_state()


# ── Helpers ────────────────────────────────────────────────────────────────────

def _variant_dataframe(variants):
    import pandas as pd
    rows = []
    for v in variants:
        d = v.to_dict() if hasattr(v, "to_dict") else v
        rows.append({
            "Gene":     d.get("gene", ""),
            "Variant":  d.get("hgvs_p") or d.get("hgvs_c") or f"{d.get('ref','')}>{d.get('alt','')}",
            "Zygosity": d.get("zygosity", ""),
            "Impact":   d.get("impact", ""),
            "ClinVar":  d.get("clinvar_significance") or "—",
            "gnomAD AF": f"{d['gnomad_af']:.5f}" if d.get("gnomad_af") is not None else "—",
            "Position": f"{d.get('chrom','')}:{d.get('pos','')}",
            "rsID":     d.get("rsid", "—") if d.get("rsid") != "." else "—",
        })
    return pd.DataFrame(rows)


def _get_parser() -> Optional[VCFParser]:
    path = st.session_state.get("vcf_path")
    if path and Path(path).exists():
        try:
            return VCFParser(path)
        except Exception:
            pass
    return None


def _get_analyzer() -> Optional[GenomeAnalyzer]:
    key = get_api_key(username=_AUTH_USER) or os.environ.get("ANTHROPIC_API_KEY", "")
    if key:
        return GenomeAnalyzer(api_key=key)
    return None


def _save_current_profile() -> None:
    """Persist session panel results + chat to the current profile."""
    name = st.session_state.get("current_profile")
    if not name:
        return
    profile = load_profile(name, username=_AUTH_USER)
    # Merge panel results
    for panel_name, analysis in st.session_state.get("panel_analyses", {}).items():
        variants = st.session_state.get("panel_variants", {}).get(panel_name, [])
        profile["panel_results"][panel_name] = {
            "variants": serialize_panel_variants(variants),
            "analysis": analysis,
            "date": date.today().isoformat(),
        }
    # Save chat
    profile["chat_api_messages"] = serialize_api_messages(
        st.session_state.get("chat_api_messages", [])
    )
    profile["chat_display"] = serialize_chat_display(
        st.session_state.get("chat_display", [])
    )
    profile["vcf_path"] = st.session_state.get("vcf_path")
    profile["ancestry_analysis"] = st.session_state.get("ancestry_analysis", "")
    profile["ancestry_rsid_variants"] = serialize_panel_variants(
        st.session_state.get("ancestry_rsid_variants", [])
    )
    profile["ancestry_gene_variants"] = serialize_panel_variants(
        st.session_state.get("ancestry_gene_variants", [])
    )
    save_profile(name, profile, username=_AUTH_USER)


def _load_profile_into_state(name: str) -> None:
    """Load saved profile data into session state."""
    profile = load_profile(name, username=_AUTH_USER)
    st.session_state["current_profile"] = name
    st.session_state["vcf_path"] = profile.get("vcf_path")

    # Restore panel results
    panel_variants: dict = {}
    panel_analyses: dict = {}
    panel_results = profile.get("panel_results", {})
    for panel_name, result in panel_results.items():
        panel_variants[panel_name] = result.get("variants", [])
        panel_analyses[panel_name] = result.get("analysis", "")

    st.session_state["panel_variants"] = panel_variants
    st.session_state["panel_analyses"] = panel_analyses
    st.session_state["analysis_done"] = bool(panel_analyses)

    # Restore chat
    st.session_state["chat_api_messages"] = profile.get("chat_api_messages", [])
    st.session_state["chat_display"] = profile.get("chat_display", [])
    st.session_state["summary_text"] = ""
    st.session_state["ancestry_analysis"] = profile.get("ancestry_analysis", "")
    st.session_state["ancestry_rsid_variants"] = profile.get("ancestry_rsid_variants", [])
    st.session_state["ancestry_gene_variants"] = profile.get("ancestry_gene_variants", [])


def _stream_to_placeholder(gen, placeholder) -> str:
    full = ""
    for chunk in gen:
        full += chunk
        placeholder.markdown(full + "▌")
    placeholder.markdown(full)
    return full


# ── Sidebar ────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.title("🧬 Genomic VCF Analyzer")
    st.caption("Whole-genome variant analysis · Claude AI")

    # ── Logged-in user + logout ───────────────────────────────────────────────
    col_user, col_logout = st.columns([3, 1])
    with col_user:
        st.caption(f"Logged in as **{_AUTH_USER}**")
    with col_logout:
        if st.button("↩", help="Log out", use_container_width=True):
            for key in list(st.session_state.keys()):
                del st.session_state[key]
            st.rerun()

    st.divider()

    # ── 1. Profile management ──────────────────────────────────────────────────
    st.subheader("👤 Profile")
    profiles = list_profiles(username=_AUTH_USER)

    col_sel, col_new = st.columns([3, 1])
    with col_sel:
        profile_options = profiles if profiles else []
        current_idx = (
            profile_options.index(st.session_state["current_profile"])
            if st.session_state["current_profile"] in profile_options
            else 0
        )
        selected_profile = st.selectbox(
            "Select profile",
            options=profile_options,
            index=current_idx if profile_options else 0,
            placeholder="No profiles yet",
            label_visibility="collapsed",
        ) if profile_options else None

    with col_new:
        new_profile_btn = st.button("＋", help="Create new profile", use_container_width=True)

    if new_profile_btn:
        st.session_state["creating_profile"] = True

    if st.session_state.get("creating_profile"):
        new_name = st.text_input("New profile name", key="new_profile_name_input",
                                  placeholder="e.g. John Smith")
        c1, c2 = st.columns(2)
        with c1:
            if st.button("Create", use_container_width=True) and new_name.strip():
                name = new_name.strip()
                # Save to disk first so list_profiles() includes it on rerun
                save_profile(name, load_profile(name, username=_AUTH_USER), username=_AUTH_USER)
                _load_profile_into_state(name)
                st.session_state["creating_profile"] = False
                st.session_state["show_vcf_upload"] = True
                st.rerun()
        with c2:
            if st.button("Cancel", use_container_width=True):
                st.session_state["creating_profile"] = False
                st.rerun()

    if selected_profile and selected_profile != st.session_state.get("current_profile"):
        _load_profile_into_state(selected_profile)
        st.rerun()

    if selected_profile:
        col_del, col_save = st.columns(2)
        with col_del:
            if st.button("🗑 Delete", use_container_width=True,
                         help="Delete this profile and its data"):
                delete_profile(selected_profile, username=_AUTH_USER)
                st.session_state["current_profile"] = None
                st.session_state["vcf_path"] = None
                st.rerun()
        with col_save:
            if st.button("💾 Save", use_container_width=True,
                         help="Save current session to profile"):
                _save_current_profile()
                st.success("Saved!")

    st.divider()

    # ── 2. Genome file ─────────────────────────────────────────────────────────
    st.subheader("📂 Genome File")

    vcf_path = st.session_state.get("vcf_path")
    vcf_exists = vcf_path and Path(vcf_path).exists()

    if vcf_exists and not st.session_state.get("show_vcf_upload"):
        st.success(f"✓ {Path(vcf_path).name}")
        has_tbi = Path(vcf_path + ".tbi").exists()
        if has_tbi:
            st.caption("✓ Tabix index present")
        else:
            st.caption("⚠️ No tabix index")
            if tabix_available():
                if st.button("⚡ Create index now"):
                    with st.spinner("Indexing VCF (may take a few minutes)…"):
                        ok = create_tabix_index(vcf_path)
                    if ok:
                        st.success("Index created!")
                        st.rerun()
                    else:
                        st.error("Indexing failed.")

        # Show detected genome build and VCF sanity diagnostics
        try:
            _diag_parser = VCFParser(vcf_path)
            diag = _diag_parser.get_diagnostics()
            build_label = diag["genome_build"].upper()
            build_color = "🟢" if diag["genome_build"] == "hg38" else "🟡"
            st.caption(f"{build_color} Build: **{build_label}** · {diag['chrom_style']}")
            if diag["genome_build"] == "hg19":
                st.warning(
                    "hg19/GRCh37 detected — coordinates will be fetched for this build.",
                )
            ann_fmt = diag["annotation_format"].upper() if diag["annotation_format"] != "none" else "none"
            st.caption(f"Annotations: {ann_fmt} · Samples: {diag['num_samples']}")

            # Chromosome names from the tabix index
            if diag["seqnames_sample"]:
                names_str = ", ".join(diag["seqnames_sample"])
                if len(_diag_parser._seqnames) > 8:
                    names_str += f", … ({len(_diag_parser._seqnames)} total)"
                st.caption(f"Chromosomes: `{names_str}`")

            # Sample variants — confirms the file is actually readable
            if diag["sample_variants"]:
                st.caption("First variants: " + " · ".join(diag["sample_variants"]))
            else:
                st.warning(
                    "⚠️ Could not read any variants from this file. "
                    "It may be empty, corrupted, or the wrong format.",
                )
        except Exception as _e:
            st.caption(f"⚠️ Could not read VCF diagnostics: {_e}")

        if st.button("📂 Change VCF file"):
            st.session_state["show_vcf_upload"] = True
            st.rerun()
    else:
        # ── Path input (recommended for large WGS files) ───────────────────
        st.caption("**Option A — Paste file path** (recommended for large WGS files)")
        path_input = st.text_input(
            "Full path to VCF or VCF.gz",
            placeholder="/Users/you/data/genome.vcf.gz",
            label_visibility="collapsed",
            key="vcf_path_input",
        )
        if st.button("Use this path", disabled=not path_input.strip()):
            p = Path(path_input.strip())
            if not p.exists():
                st.error(f"File not found: {p}")
            elif not str(p).endswith((".vcf", ".vcf.gz", ".gz")):
                st.error("File must be a .vcf or .vcf.gz file.")
            else:
                use_path = str(p)

                # If the file is on a read-only location (USB, external drive)
                # and needs indexing, copy it to the profile directory first
                if str(p).endswith(".gz") and not Path(str(p) + ".tbi").exists():
                    parent_writable = os.access(str(p.parent), os.W_OK)
                    if not parent_writable:
                        profile_name = st.session_state.get("current_profile", "default")
                        dest_dir = user_profiles_dir(_AUTH_USER) / profile_name
                        dest_dir.mkdir(parents=True, exist_ok=True)
                        dest = dest_dir / p.name
                        if not dest.exists():
                            import shutil
                            with st.spinner(f"Copying VCF to local storage ({p.stat().st_size // (1024*1024)} MB)…"):
                                shutil.copy2(str(p), str(dest))
                            st.info(f"Copied to {dest} (original location is read-only)")
                        use_path = str(dest)

                    if tabix_available():
                        with st.spinner("Creating tabix index (this may take a few minutes)…"):
                            ok = create_tabix_index(use_path)
                        if ok:
                            st.success("✓ Index created!")
                        else:
                            st.warning("⚠️ Could not create tabix index. Queries will be slow.")
                    else:
                        st.warning("⚠️ No tabix index found. Install htslib: `brew install htslib`")

                st.session_state["vcf_path"] = use_path
                st.session_state["show_vcf_upload"] = False
                _save_current_profile()
                st.rerun()

        st.divider()

        # ── File uploader (for smaller files only) ─────────────────────────
        st.caption("**Option B — Upload file** (up to 1 GB)")
        vcf_upload = st.file_uploader(
            "VCF or VCF.gz file",
            type=["gz", "vcf"],
            help="Upload limit is 1 GB. For very large files or files on external drives, use Option A above.",
            label_visibility="collapsed",
        )
        if vcf_upload:
            profile_name = st.session_state.get("current_profile", "default")
            with st.spinner(f"Saving {vcf_upload.name}…"):
                vcf_bytes = vcf_upload.read()
                saved_path, tbi_ok = save_vcf_for_profile(
                    profile_name, vcf_bytes, vcf_upload.name, username=_AUTH_USER
                )
                st.session_state["vcf_path"] = saved_path

            if tbi_ok:
                st.success("✓ Saved and indexed!")
            elif vcf_upload.name.endswith(".gz"):
                st.warning("⚠️ Saved but could not create tabix index. Install htslib: `brew install htslib`")
            else:
                st.success("✓ Saved!")

            st.session_state["show_vcf_upload"] = False
            _save_current_profile()
            st.rerun()

    st.divider()

    # ── 3. API key ─────────────────────────────────────────────────────────────
    st.subheader("🔑 API Key")
    saved_key = get_api_key(username=_AUTH_USER)
    api_key = saved_key or os.environ.get("ANTHROPIC_API_KEY", "")

    if st.session_state["editing_api_key"] or not api_key:
        new_key = st.text_input(
            "Anthropic API Key",
            type="password",
            placeholder="sk-ant-…",
            label_visibility="collapsed",
        )
        col_save_key, col_cancel_key = st.columns(2)
        with col_save_key:
            if st.button("Save key", use_container_width=True) and new_key.strip():
                save_api_key(new_key.strip(), username=_AUTH_USER)
                api_key = new_key.strip()
                st.session_state["editing_api_key"] = False
                st.rerun()
        with col_cancel_key:
            if api_key and st.button("Cancel", use_container_width=True):
                st.session_state["editing_api_key"] = False
                st.rerun()
    else:
        masked = api_key[:10] + "••••••••" + api_key[-3:] if len(api_key) > 13 else "••••••••••"
        st.code(masked, language=None)
        if st.button("✏ Change key", use_container_width=True):
            st.session_state["editing_api_key"] = True
            st.rerun()

    st.divider()

    # ── 4. Panel selection + run ───────────────────────────────────────────────
    st.subheader("🔬 Select Panels")
    selected_panels: list[str] = []
    for panel_name, gene_set in ALL_PANELS.items():
        icon = PANEL_ICONS.get(panel_name, "🔹")
        checked = st.checkbox(
            f"{icon} {panel_name.replace('_', ' ').title()} ({len(gene_set)})",
            value=(panel_name in ("cancer", "cardiovascular")),
            help=PANEL_DESCRIPTIONS.get(panel_name, ""),
        )
        if checked:
            selected_panels.append(panel_name)

    st.divider()
    run_disabled = not vcf_exists or not api_key or not selected_panels
    run_btn = st.button("▶ Run Panel Analysis", type="primary",
                        use_container_width=True, disabled=run_disabled)
    if not vcf_exists:
        st.caption("⬆ Upload a VCF file above.")
    elif not api_key:
        st.caption("🔑 Enter API key above.")
    elif not selected_panels:
        st.caption("☑ Select at least one panel.")

    # ── Admin panel ───────────────────────────────────────────────────────────
    st.divider()
    with st.expander("🔧 Admin", expanded=False):
        if st.session_state.get("admin_verified"):
            st.success("Admin access granted")
            users = list_users()
            st.caption(f"**{len(users)} registered user(s)**")
            for u in users:
                col_name, col_del = st.columns([3, 1])
                with col_name:
                    st.text(f"{u['username']} (joined {u['created'][:10]})")
                with col_del:
                    if st.button("🗑", key=f"del_user_{u['username']}",
                                 help=f"Delete {u['username']}"):
                        if u["username"] == _AUTH_USER:
                            st.error("Cannot delete yourself.")
                        else:
                            delete_user(u["username"])
                            st.success(f"Deleted {u['username']}")
                            st.rerun()
            st.divider()
            st.caption("**Reset a user's password**")
            reset_target = st.selectbox(
                "User", [u["username"] for u in users],
                key="reset_target", label_visibility="collapsed",
            )
            new_pw = st.text_input("New password", type="password", key="reset_pw")
            if st.button("Reset password") and reset_target and new_pw:
                ok, msg = reset_user_password(reset_target, new_pw)
                if ok:
                    st.success(msg)
                else:
                    st.error(msg)
        else:
            admin_pw = st.text_input("Admin password", type="password", key="admin_pw_input")
            if st.button("Unlock"):
                if verify_admin(admin_pw):
                    st.session_state["admin_verified"] = True
                    st.rerun()
                else:
                    st.error("Incorrect admin password.")


# ── Disclaimer ─────────────────────────────────────────────────────────────────
st.markdown("""
<div class="disclaimer-box">
⚠️ <strong>Disclaimer</strong> — Research and educational purposes only.
Not medical advice. Consult qualified healthcare professionals for clinical decisions.
Significant variants should be confirmed by a CLIA-certified laboratory.
</div>
""", unsafe_allow_html=True)

# Profile header
if st.session_state.get("current_profile"):
    st.markdown(
        f'<span class="profile-badge">👤 {st.session_state["current_profile"]}</span>',
        unsafe_allow_html=True,
    )
    st.markdown("")

# ── Tabs ───────────────────────────────────────────────────────────────────────
tab_analysis, tab_query, tab_chat, tab_ancestry, tab_summary, tab_about = st.tabs([
    "📊 Panel Analysis", "💬 Quick Query", "🔬 Genomic Chat",
    "🌍 Ancestry", "📋 Summary & Export", "ℹ️ About",
])


# ══ Tab 1: Panel Analysis ══════════════════════════════════════════════════════
with tab_analysis:

    if run_btn and vcf_exists and api_key and selected_panels:
        try:
            parser = VCFParser(vcf_path)
            analyzer = GenomeAnalyzer(api_key=api_key)
        except Exception as e:
            st.error(f"Initialization error: {e}")
            st.stop()

        if not parser._has_tabix:
            st.warning("⚠️ No tabix index. Analysis may be slow for WGS files.")

        for panel_name in selected_panels:
            gene_set = ALL_PANELS[panel_name]
            icon = PANEL_ICONS.get(panel_name, "🔹")
            st.markdown(f"## {icon} {panel_name.replace('_',' ').title()} Panel")

            with st.spinner(f"Scanning {len(gene_set)} genes…"):
                try:
                    variants = parser.filter_panel_variants(panel_name, gene_set)
                except Exception as e:
                    st.error(f"Error: {e}")
                    st.session_state["panel_variants"][panel_name] = []
                    st.session_state["panel_analyses"][panel_name] = ""
                    continue

            st.session_state["panel_variants"][panel_name] = variants

            if variants:
                st.success(f"Found **{len(variants)}** clinically filtered variant(s)")
                st.dataframe(_variant_dataframe(variants), use_container_width=True, hide_index=True)
            else:
                st.info(f"No clinically relevant variants found — Claude will provide context.")

            st.markdown("#### 🤖 Claude's Interpretation")
            ph = st.empty()
            with st.spinner("Analyzing…"):
                try:
                    full = _stream_to_placeholder(
                        analyzer.analyze_panel(variants, panel_name, gene_set), ph
                    )
                    st.session_state["panel_analyses"][panel_name] = full
                except Exception as e:
                    st.error(f"Claude API error: {e}")
                    st.session_state["panel_analyses"][panel_name] = ""
            st.divider()

        st.session_state["analysis_done"] = True
        _save_current_profile()

        # Download full report
        if st.session_state["panel_analyses"]:
            import tempfile
            with tempfile.NamedTemporaryFile(delete=False, suffix=".md", mode="w") as tmp:
                rpath = tmp.name
            save_markdown_report(
                vcf_path=vcf_path or "",
                panels_analyzed=selected_panels,
                panel_variants=st.session_state["panel_variants"],
                panel_analyses=st.session_state["panel_analyses"],
                output_path=rpath,
            )
            st.download_button(
                "⬇ Download Full Report (.md)",
                data=Path(rpath).read_bytes(),
                file_name=f"genomic_report_{date.today()}.md",
                mime="text/markdown",
                type="primary",
            )

    elif st.session_state["analysis_done"] and st.session_state["panel_analyses"]:
        # Restore saved analysis
        for panel_name, analysis in st.session_state["panel_analyses"].items():
            icon = PANEL_ICONS.get(panel_name, "🔹")
            st.markdown(f"## {icon} {panel_name.replace('_',' ').title()} Panel")
            variants = st.session_state["panel_variants"].get(panel_name, [])
            if variants:
                st.dataframe(_variant_dataframe(variants), use_container_width=True, hide_index=True)
            if analysis:
                st.markdown("#### 🤖 Claude's Interpretation")
                st.markdown(analysis)
            st.divider()
    else:
        st.markdown("""
        ### How to run a panel analysis
        1. Select or create a **profile** in the sidebar
        2. Upload your **VCF or VCF.gz** file (tabix index created automatically)
        3. Enter your **Anthropic API key** (saved for future sessions)
        4. Check the **panels** you want to analyze
        5. Click **▶ Run Panel Analysis**

        Previously saved results for this profile will appear here automatically.
        """)


# ══ Tab 2: Quick Query ════════════════════════════════════════════════════════
with tab_query:
    st.markdown("### Ask a specific question about this genome")
    st.caption("For multi-turn conversation, use the **🔬 Genomic Chat** tab instead.")

    question = st.text_input("Your question", placeholder="What is my APOE genotype?",
                              label_visibility="collapsed")
    ask_btn = st.button("🔍 Ask Claude", type="primary",
                         disabled=(not question or not vcf_exists or not api_key))

    if not vcf_exists:
        st.caption("⬆ Upload a VCF file in the sidebar.")
    elif not api_key:
        st.caption("🔑 Enter API key in the sidebar.")

    if ask_btn and question and vcf_exists and api_key:
        try:
            analyzer = GenomeAnalyzer(api_key=api_key)
            parser = VCFParser(vcf_path)
        except Exception as e:
            st.error(f"Initialization error: {e}")
            st.stop()

        with st.spinner("Identifying relevant genes…"):
            try:
                genes = analyzer.identify_relevant_genes(question)
            except Exception as e:
                st.error(f"Gene identification error: {e}")
                st.stop()

        from gene_lookup import get_coordinates_batch
        coord_map = get_coordinates_batch(genes, build=parser.genome_build)
        found_genes = [g for g, c in coord_map.items() if c is not None]
        not_found = [g for g, c in coord_map.items() if c is None]

        if not found_genes:
            st.error("Could not locate any of the identified genes. Try the Genomic Chat tab for more flexible queries.")
            st.stop()

        st.info(f"Querying: **{', '.join(found_genes)}**")
        if not_found:
            st.caption(f"Could not resolve: {', '.join(not_found)}")

        with st.spinner("Retrieving variants…"):
            try:
                variants = parser.get_variants_for_question(found_genes)
            except Exception as e:
                st.error(f"VCF error: {e}")
                st.stop()

        st.session_state["query_variants"] = variants
        if variants:
            st.dataframe(_variant_dataframe(variants), use_container_width=True, hide_index=True)
        else:
            st.info("No clinically filtered variants found — reference genotype likely.")

        st.markdown("#### 🤖 Claude's Answer")
        ph = st.empty()
        try:
            full = _stream_to_placeholder(
                analyzer.answer_question(variants, question, context_genes=found_genes), ph
            )
            st.session_state["query_answer"] = full
        except Exception as e:
            st.error(f"Claude API error: {e}")

        if st.session_state.get("query_answer"):
            st.download_button(
                "⬇ Download Answer",
                data=f"# Genomic Query\n\n**Question**: {question}\n\n{st.session_state['query_answer']}\n\n*Research use only*",
                file_name=f"query_{date.today()}.md",
                mime="text/markdown",
            )

    elif st.session_state.get("query_answer"):
        variants = st.session_state.get("query_variants", [])
        if variants:
            st.dataframe(_variant_dataframe(variants), use_container_width=True, hide_index=True)
        st.markdown("#### 🤖 Claude's Answer")
        st.markdown(st.session_state["query_answer"])


# ══ Tab 3: Genomic Chat ════════════════════════════════════════════════════════
with tab_chat:
    st.markdown("### 🔬 Genomic Chat")
    st.markdown(
        "Ask anything in plain English — Claude identifies relevant genes, queries your VCF, "
        "and interprets findings. The conversation persists across sessions for this profile."
    )
    st.caption(
        "Examples: *What are my estrogen receptor variants?* · "
        "*Check my asthma genes* · "
        "*Do I have pharmacogenomic variants affecting antidepressants?*"
    )

    chat_ready = vcf_exists and bool(api_key)

    # ── Scrollable messages area ───────────────────────────────────────────────
    # Defined first so it renders above the input bar.  Filled below.
    chat_area = st.container(height=700)

    # ── Input bar — always visible, fixed below the messages area ─────────────
    user_input = st.chat_input("Ask about this genome…", disabled=not chat_ready)

    if not chat_ready:
        if not vcf_exists:
            st.info("⬆ Upload a VCF file in the sidebar to start chatting.")
        elif not api_key:
            st.info("🔑 Enter API key in the sidebar to start.")

    # ── All message rendering lives inside the scrollable container ───────────
    with chat_area:

        # Render existing chat history
        for entry in st.session_state["chat_display"]:
            with st.chat_message(entry["role"]):
                # Database lookup events (stored under "db_lookups")
                for dl in entry.get("db_lookups", []):
                    dbs = dl.get("databases", [])
                    dgenes = dl.get("genes", [])
                    rsids_found = dl.get("rsids_found", [])
                    counts = dl.get("counts", {})
                    db_label = ", ".join(s.replace("_", " ").title() for s in dbs)
                    with st.expander(
                        f"🗄️ Live DB lookup — {db_label}: {', '.join(dgenes)}",
                        expanded=False,
                    ):
                        st.caption(dl.get("reason", ""))
                        for g, srcs in counts.items():
                            parts = [f"{s}: {n}" for s, n in srcs.items()]
                            st.caption(f"**{g}** — " + " | ".join(parts))
                        if rsids_found:
                            st.caption(f"🔑 rsIDs discovered: {', '.join(rsids_found[:20])}"
                                       + (" …" if len(rsids_found) > 20 else ""))
                        else:
                            st.caption("No rsIDs found in databases for these genes.")
                # VCF tool call events (stored under "tool_calls")
                for tc in entry.get("tool_calls", []):
                    genes = tc.get("genes", [])
                    variants = tc.get("variants", [])
                    fetched = tc.get("fetched_online", [])
                    not_in_db = tc.get("not_in_db", [])
                    with st.expander(f"🔍 Queried {len(genes)} gene(s): {', '.join(genes)}",
                                      expanded=False):
                        st.caption(tc.get("reason", ""))
                        if fetched:
                            st.caption(f"🌐 Coordinates fetched from mygene.info: {', '.join(fetched)}")
                        if not_in_db:
                            st.caption(f"⚠️ No coordinates found: {', '.join(not_in_db)}")
                        if variants:
                            rows = [v.to_dict() if hasattr(v, "to_dict") else v for v in variants]
                            if rows:
                                import pandas as pd
                                st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
                        else:
                            st.caption("No clinically filtered variants — reference genotype.")
                if entry.get("text"):
                    st.markdown(entry["text"])
                if entry.get("error"):
                    st.error(entry["error"])

        # ── Live exchange for the current turn ────────────────────────────────
        if user_input and chat_ready:
            try:
                import anthropic as _anthropic
                client = _anthropic.Anthropic(api_key=api_key)
                parser = VCFParser(vcf_path)
            except Exception as e:
                st.error(f"Initialization error: {e}")
                st.stop()

            st.session_state["chat_api_messages"].append(
                {"role": "user", "content": user_input}
            )

            with st.chat_message("user"):
                st.markdown(user_input)

            display_entry: dict = {
                "role": "assistant",
                "db_lookups": [],
                "tool_calls": [],
                "text": "",
                "error": "",
            }

            with st.chat_message("assistant"):
                text_ph = st.empty()
                accumulated = ""

                for event in run_chat_turn(
                    client=client, parser=parser,
                    messages=st.session_state["chat_api_messages"],
                ):
                    if isinstance(event, DatabaseLookupEvent):
                        db_label = ", ".join(s.replace("_"," ").title()
                                             for s in event.databases)
                        display_entry["db_lookups"].append({
                            "genes": event.genes, "databases": event.databases,
                            "reason": event.reason, "counts": {}, "rsids_found": [],
                        })
                        with st.expander(
                            f"🗄️ Querying live databases — {db_label}: {', '.join(event.genes)}",
                            expanded=True,
                        ):
                            st.caption(event.reason)
                            st.caption("⏳ Fetching from ClinVar / GWAS Catalog / PharmGKB…")

                    elif isinstance(event, DatabaseResultEvent):
                        if display_entry["db_lookups"]:
                            dl = display_entry["db_lookups"][-1]
                            dl["counts"] = event.counts
                            dl["rsids_found"] = event.rsids_found
                        db_label = ", ".join(s.replace("_"," ").title()
                                             for s in event.databases)
                        with st.expander(
                            f"🗄️ DB results — {db_label}: {', '.join(event.genes)}",
                            expanded=True,
                        ):
                            for gene, srcs in event.counts.items():
                                parts = [f"{s}: {n}" for s, n in srcs.items()]
                                st.caption(f"**{gene}** — " + " | ".join(parts))
                            if event.rsids_found:
                                st.caption(
                                    f"🔑 {len(event.rsids_found)} rsIDs discovered: "
                                    + ", ".join(event.rsids_found[:20])
                                    + (" …" if len(event.rsids_found) > 20 else "")
                                )
                            else:
                                st.caption("No rsIDs found — will query by gene region.")

                    elif isinstance(event, ToolCallEvent):
                        display_entry["tool_calls"].append({
                            "genes": event.genes, "reason": event.reason,
                            "variants": [], "fetched_online": event.fetched_online, "not_in_db": [],
                        })
                        online_note = (
                            f" · 🌐 fetching {len(event.fetched_online)} via mygene.info…"
                            if event.fetched_online else ""
                        )
                        with st.expander(
                            f"🔍 Querying {len(event.genes)} gene(s): {', '.join(event.genes)}",
                            expanded=True,
                        ):
                            st.caption(event.reason)
                            st.caption(f"⏳ Querying VCF…{online_note}")

                    elif isinstance(event, ToolResultEvent):
                        if display_entry["tool_calls"]:
                            tc = display_entry["tool_calls"][-1]
                            tc["variants"] = event.variants
                            tc["not_in_db"] = event.not_in_db
                            tc["fetched_online"] = event.fetched_online
                        with st.expander(
                            f"🔍 Found {len(event.variants)} variant(s) in: {', '.join(event.genes)}",
                            expanded=len(event.variants) > 0,
                        ):
                            if event.fetched_online:
                                st.caption(f"🌐 mygene.info: {', '.join(event.fetched_online)}")
                            if event.not_in_db:
                                st.caption(f"⚠️ No coords found: {', '.join(event.not_in_db)}")
                            if event.variants:
                                st.dataframe(
                                    _variant_dataframe(event.variants),
                                    use_container_width=True, hide_index=True,
                                )
                            else:
                                st.caption("No clinically filtered variants — reference genotype.")

                    elif isinstance(event, TextEvent):
                        accumulated += event.text
                        text_ph.markdown(accumulated)
                        display_entry["text"] = accumulated

                    elif isinstance(event, ErrorEvent):
                        st.error(event.message)
                        display_entry["error"] = event.message

            st.session_state["chat_display"].append(
                {"role": "user", "text": user_input, "tool_calls": []}
            )
            st.session_state["chat_display"].append(display_entry)
            _save_current_profile()

    # ── Controls below the chat area ──────────────────────────────────────────
    if st.session_state["chat_display"]:
        if st.button("🗑 Clear conversation", key="clear_chat"):
            st.session_state["chat_api_messages"] = []
            st.session_state["chat_display"] = []
            _save_current_profile()
            st.rerun()


# ══ Tab 4: Ancestry ═══════════════════════════════════════════════════════════
with tab_ancestry:
    st.markdown("### 🌍 Genomic Ancestry Analysis")
    st.markdown(
        "Queries **ancestry-informative markers (AIMs)** and the **full mitochondrial genome** "
        "to give directional signals about genomic ancestry, maternal haplogroup, and "
        "related phenotypic traits."
    )
    st.info(
        "**What this covers:** "
        "(1) **mtDNA haplogroup** — maternal lineage (H, J, K, L, M, …) from a full chrM scan; "
        "(2) **Continental ancestry signals** — ~20 autosomal AIMs (pigmentation, diet, "
        "morphology, immune adaptation) plus surrounding gene regions. "
        "Provides *directional* signals, not percentage breakdowns. "
        "For precise admixture proportions a dedicated tool (e.g. ADMIXTURE with a "
        "1000 Genomes reference panel) is required. Genomic ancestry ≠ ethnic identity.",
        icon="ℹ️",
    )

    ancestry_ready = vcf_exists and bool(api_key)
    if not vcf_exists:
        st.info("⬆ Upload a VCF file in the sidebar to run ancestry analysis.")
    elif not api_key:
        st.info("🔑 Enter API key in the sidebar.")

    if ancestry_ready:
        run_ancestry_btn = st.button(
            "🌍 Run Ancestry Analysis",
            type="primary",
            disabled=not ancestry_ready,
        )

        if run_ancestry_btn:
            try:
                parser = VCFParser(vcf_path)
                analyzer = GenomeAnalyzer(api_key=api_key)
            except Exception as e:
                st.error(f"Initialization error: {e}")
                st.stop()

            # ── Step 1: rsID direct lookup for key AIMs ────────────────────
            st.markdown("#### Step 1 — Querying known ancestry-informative rsIDs")
            with st.spinner(f"Looking up {len(ANCESTRY_RSIDS)} ancestry-informative SNPs…"):
                try:
                    rsid_variants = parser.get_variants_by_rsids(ANCESTRY_RSIDS)
                except Exception as e:
                    st.error(f"rsID lookup error: {e}")
                    rsid_variants = []

            st.success(f"Found **{len(rsid_variants)}** variants at ancestry rsID positions")
            if rsid_variants:
                st.dataframe(
                    _variant_dataframe(rsid_variants),
                    use_container_width=True, hide_index=True,
                )

            # ── Step 2: unfiltered gene-region query ───────────────────────
            st.markdown("#### Step 2 — Scanning ancestry gene regions (unfiltered)")
            with st.spinner(f"Scanning {len(ANCESTRY_GENES)} ancestry-informative genes…"):
                try:
                    gene_variants = parser.get_unfiltered_variants(
                        list(ANCESTRY_GENES), max_per_gene=100
                    )
                except Exception as e:
                    st.error(f"Gene region query error: {e}")
                    gene_variants = []

            # Deduplicate against rsid_variants
            seen_keys = {(v.chrom, v.pos, v.ref, v.alt) for v in rsid_variants}
            gene_variants = [
                v for v in gene_variants
                if (v.chrom, v.pos, v.ref, v.alt) not in seen_keys
            ]
            st.success(f"Found **{len(gene_variants)}** additional variants in ancestry gene regions")
            if gene_variants:
                st.dataframe(
                    _variant_dataframe(gene_variants[:50]),
                    use_container_width=True, hide_index=True,
                )
                if len(gene_variants) > 50:
                    st.caption(f"Showing top 50 of {len(gene_variants)} gene-region variants.")

            # ── Step 3: Claude interpretation ─────────────────────────────
            st.markdown("#### 🤖 Claude's Ancestry Interpretation")
            ph = st.empty()
            with st.spinner("Analysing ancestry markers…"):
                try:
                    full = _stream_to_placeholder(
                        analyzer.analyze_ancestry(rsid_variants, gene_variants), ph
                    )
                    st.session_state["ancestry_analysis"] = full
                    st.session_state["ancestry_rsid_variants"] = rsid_variants
                    st.session_state["ancestry_gene_variants"] = gene_variants
                    _save_current_profile()
                except Exception as e:
                    st.error(f"Claude API error: {e}")

        elif st.session_state.get("ancestry_analysis"):
            # Show previously computed ancestry analysis
            rsid_v = st.session_state.get("ancestry_rsid_variants", [])
            gene_v = st.session_state.get("ancestry_gene_variants", [])
            if rsid_v:
                st.markdown("##### Ancestry-informative SNPs found")
                st.dataframe(_variant_dataframe(rsid_v), use_container_width=True, hide_index=True)
            if gene_v:
                st.markdown("##### Additional gene-region variants")
                st.dataframe(
                    _variant_dataframe(gene_v[:50]),
                    use_container_width=True, hide_index=True,
                )
            st.markdown("#### 🤖 Claude's Ancestry Interpretation")
            st.markdown(st.session_state["ancestry_analysis"])

            # Download button
            profile_nm = st.session_state.get("current_profile", "patient")
            st.download_button(
                "⬇ Download Ancestry Report",
                data=(
                    f"# Genomic Ancestry Analysis — {profile_nm}\n"
                    f"*Generated: {date.today().isoformat()}*\n\n"
                    f"{st.session_state['ancestry_analysis']}\n\n"
                    "---\n*Genomic ancestry analysis. For research and educational "
                    "purposes only. Not medical advice.*"
                ),
                file_name=f"ancestry_{profile_nm}_{date.today()}.md",
                mime="text/markdown",
            )


# ══ Tab 5: Summary & Export ════════════════════════════════════════════════════
with tab_summary:
    profile_name = st.session_state.get("current_profile", "Patient")
    st.markdown(f"### 📋 Summary — {profile_name}")

    has_panel = bool(st.session_state.get("panel_analyses"))
    has_chat = bool(st.session_state.get("chat_display"))

    if not has_panel and not has_chat:
        st.info(
            "No results yet for this profile. Run a **Panel Analysis** or use "
            "**Genomic Chat** first, then come back here to generate a summary."
        )
    else:
        # Show what's available to summarize
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Panels analyzed", len(st.session_state.get("panel_analyses", {})))
        with col2:
            chat_turns = sum(
                1 for e in st.session_state.get("chat_display", [])
                if e.get("role") == "user"
            )
            st.metric("Chat queries", chat_turns)

        st.divider()

        if st.session_state.get("summary_text"):
            st.markdown(st.session_state["summary_text"])
        else:
            st.markdown(
                "Click **Generate Summary** to have Claude synthesize all findings "
                "from your panel analyses and chat queries into a single clinical report."
            )

        col_gen, col_dl = st.columns([2, 1])
        with col_gen:
            gen_btn = st.button(
                "✨ Generate AI Summary",
                type="primary",
                disabled=not api_key,
                help="Uses Claude to synthesize all findings into a structured report.",
            )
        with col_dl:
            if st.session_state.get("summary_text"):
                profile_nm = st.session_state.get("current_profile", "patient")
                summary_content = (
                    f"# Genomic Summary — {profile_nm}\n"
                    f"*Generated: {date.today().isoformat()}*\n\n"
                    f"{st.session_state['summary_text']}\n\n"
                    f"---\n*Research and educational purposes only. Not medical advice.*"
                )
                st.download_button(
                    "⬇ Download Summary",
                    data=summary_content,
                    file_name=f"genomic_summary_{profile_nm}_{date.today()}.md",
                    mime="text/markdown",
                    use_container_width=True,
                )

        if gen_btn and api_key:
            try:
                analyzer = GenomeAnalyzer(api_key=api_key)
            except Exception as e:
                st.error(f"Initialization error: {e}")
                st.stop()

            # Build panel_results dict for the summarizer
            panel_results_for_summary = {}
            for panel_name, analysis in st.session_state.get("panel_analyses", {}).items():
                variants = st.session_state.get("panel_variants", {}).get(panel_name, [])
                panel_results_for_summary[panel_name] = {
                    "analysis": analysis,
                    "variants": variants,
                }

            st.markdown("#### 🤖 Claude's Summary")
            ph = st.empty()
            with st.spinner("Generating summary…"):
                try:
                    full = _stream_to_placeholder(
                        analyzer.generate_summary(
                            patient_name=profile_name,
                            panel_results=panel_results_for_summary,
                            chat_display=st.session_state.get("chat_display", []),
                        ),
                        ph,
                    )
                    st.session_state["summary_text"] = full
                    _save_current_profile()
                    st.rerun()
                except Exception as e:
                    st.error(f"Claude API error: {e}")

        # Also offer full report download if panel analysis was done
        if has_panel and vcf_exists:
            st.divider()
            st.markdown("**Full detailed report** (all panel analyses with variant tables):")
            import tempfile
            with tempfile.NamedTemporaryFile(delete=False, suffix=".md", mode="w") as tmp:
                rpath = tmp.name
            try:
                save_markdown_report(
                    vcf_path=vcf_path or "",
                    panels_analyzed=list(st.session_state["panel_analyses"].keys()),
                    panel_variants=st.session_state["panel_variants"],
                    panel_analyses=st.session_state["panel_analyses"],
                    output_path=rpath,
                )
                st.download_button(
                    "⬇ Download Full Report (.md)",
                    data=Path(rpath).read_bytes(),
                    file_name=f"genomic_report_{profile_name}_{date.today()}.md",
                    mime="text/markdown",
                )
            except Exception as e:
                st.caption(f"Could not generate report: {e}")


# ══ Tab 5: About ══════════════════════════════════════════════════════════════
with tab_about:
    st.markdown(f"""
## About the Genomic VCF Analyzer

This tool analyzes whole-genome VCF files against clinical gene panels,
using **Claude claude-opus-4-7** (Anthropic) to interpret variants.

### Gene Panels

| Panel | Genes | Focus |
|-------|-------|-------|
| 🔬 Cancer | 22 | BRCA1/2, Lynch syndrome, Li-Fraumeni, APC, VHL, PTEN… |
| ❤️ Cardiovascular | 25 | Familial hypercholesterolemia, cardiomyopathies, channelopathies |
| ⏳ Longevity | 12 | APOE, FOXO3, telomere biology, sirtuins, mTOR/IGF-1 |
| ⚡ Mitochondrial | 26 | POLG, respiratory chain (nuclear); MT-ND, MT-CYB, MT-CO (mtDNA) |
| 🦴 Connective Tissue | 22 | All EDS subtypes, Loeys-Dietz, Marfan |
| 💊 Pharmacogenomics | 16 | CYP2D6/2C19/2C9, TPMT, DPYD, VKORC1, SLCO1B1, G6PD… |
| 🧠 Neurological | 22 | Parkinson's, Alzheimer's, ALS, Huntington's, ataxias |

The **Genomic Chat** tab supports queries for *any* human gene — coordinates are
fetched automatically from [mygene.info](https://mygene.info) and cached locally.

### Data storage

User profiles and settings are stored in:
`{Path.home() / ".genome_vcf_evaluator"}`

VCF files are stored per-profile in that directory. The API key is saved in
`config.json` in that directory.

### Performance tips

- **Tabix index** (`.tbi`) is created automatically on upload if `tabix` is installed.
  Without it, WGS analysis is much slower.
  Install htslib: `brew install htslib` (Mac) or `apt install tabix` (Linux).
- Gene coordinate **cache** grows over time — any gene queried via chat is cached
  for instant future lookups.

---
⚠️ **Disclaimer**: Research and educational purposes only. Not medical advice.
Consult qualified healthcare professionals for clinical decisions.
""")
