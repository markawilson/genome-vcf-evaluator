"""
Authentication for the Genome VCF Evaluator.

User credentials are stored in ~/.genome_vcf_evaluator/users.json
with PBKDF2-SHA256 hashed passwords (100k iterations + random salt).

Admin password is set via Streamlit secrets (ADMIN_PASSWORD) or
the ADMIN_PASSWORD environment variable.
"""

from __future__ import annotations

import hashlib
import json
import os
import secrets
import shutil
from datetime import datetime
from pathlib import Path

BASE_DIR = Path.home() / ".genome_vcf_evaluator"
USERS_FILE = BASE_DIR / "users.json"
PROFILES_DIR = BASE_DIR / "profiles"


def _load_users() -> dict:
    BASE_DIR.mkdir(parents=True, exist_ok=True)
    if USERS_FILE.exists():
        try:
            return json.loads(USERS_FILE.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}


def _save_users(users: dict) -> None:
    BASE_DIR.mkdir(parents=True, exist_ok=True)
    USERS_FILE.write_text(json.dumps(users, indent=2), encoding="utf-8")


def _hash_password(password: str, salt: str | None = None) -> tuple[str, str]:
    if salt is None:
        salt = secrets.token_hex(16)
    hashed = hashlib.pbkdf2_hmac(
        "sha256", password.encode("utf-8"), salt.encode("utf-8"), 100_000
    ).hex()
    return hashed, salt


def register_user(username: str, password: str) -> tuple[bool, str]:
    username = username.strip().lower()
    if not username or not password:
        return False, "Username and password are required."
    if len(username) < 3:
        return False, "Username must be at least 3 characters."
    if len(password) < 6:
        return False, "Password must be at least 6 characters."
    if not all(c.isalnum() or c in "_-" for c in username):
        return False, "Username can only contain letters, numbers, hyphens, and underscores."

    users = _load_users()
    if username in users:
        return False, "Username already taken."

    hashed, salt = _hash_password(password)
    users[username] = {
        "hash": hashed,
        "salt": salt,
        "created": datetime.now().isoformat(),
    }
    _save_users(users)
    return True, "Account created successfully."


def verify_login(username: str, password: str) -> bool:
    username = username.strip().lower()
    users = _load_users()
    if username not in users:
        return False
    user = users[username]
    hashed, _ = _hash_password(password, user["salt"])
    return hashed == user["hash"]


def get_admin_password() -> str:
    try:
        import streamlit as st
        return st.secrets["ADMIN_PASSWORD"]
    except Exception:
        return os.environ.get("ADMIN_PASSWORD", "")


def verify_admin(password: str) -> bool:
    admin_pw = get_admin_password()
    if not admin_pw:
        return False
    return password == admin_pw


def list_users() -> list[dict]:
    users = _load_users()
    return [
        {"username": u, "created": info.get("created", "unknown")}
        for u, info in sorted(users.items())
    ]


def delete_user(username: str) -> bool:
    username = username.strip().lower()
    users = _load_users()
    if username not in users:
        return False
    del users[username]
    _save_users(users)
    user_dir = PROFILES_DIR / username
    if user_dir.exists():
        shutil.rmtree(user_dir, ignore_errors=True)
    return True


def reset_user_password(username: str, new_password: str) -> tuple[bool, str]:
    username = username.strip().lower()
    if len(new_password) < 6:
        return False, "Password must be at least 6 characters."
    users = _load_users()
    if username not in users:
        return False, "User not found."
    hashed, salt = _hash_password(new_password)
    users[username]["hash"] = hashed
    users[username]["salt"] = salt
    _save_users(users)
    return True, "Password reset successfully."
