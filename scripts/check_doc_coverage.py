#!/usr/bin/env python3
"""Check MOCHA doc coverage: exports vs man aliases vs rendered pkgdown reference."""
from __future__ import annotations

import argparse
import html
import os
import re
import subprocess
import sys
from pathlib import Path


def git_show(repo: Path, ref: str, path: str) -> str:
    return subprocess.check_output(
        ["git", "show", f"{ref}:{path}"],
        cwd=repo,
        text=True,
    )


def exports_from_namespace_text(text: str) -> list[str]:
    out: list[str] = []
    for line in text.splitlines():
        m = re.match(r"^export\(([^)]+)\)", line.strip())
        if m:
            out.append(m.group(1).strip().strip("\"'"))
    return out


def rd_aliases_from_text(text: str) -> set[str]:
    return {a.strip() for a in re.findall(r"\\alias\{([^}]+)\}", text)}


def rendered_coverage(ref_dir: Path, exports: list[str]) -> set[str]:
    covered: set[str] = set()
    if not ref_dir.is_dir():
        return covered
    for path in ref_dir.glob("*.html"):
        if path.name == "index.html":
            continue
        raw = path.read_text(encoding="utf-8", errors="replace")
        txt = html.unescape(re.sub(r"<[^>]+>", " ", raw))
        txt = re.sub(r"\s+", " ", txt)
        for name in exports:
            if re.search(r"\b" + re.escape(name) + r"\s*\(", txt):
                covered.add(name)
    return covered


def collect_man_aliases(repo: Path, ref: str) -> set[str]:
    paths = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", ref, "man"],
        cwd=repo,
        text=True,
    ).splitlines()
    aliases: set[str] = set()
    for path in paths:
        if not path.endswith(".Rd"):
            continue
        aliases |= rd_aliases_from_text(git_show(repo, ref, path))
    return aliases


def collect_rendered(repo: Path, ref: str, exports: list[str]) -> set[str]:
    paths = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", ref, "docs/reference"],
        cwd=repo,
        text=True,
    ).splitlines()
    covered: set[str] = set()
    for path in paths:
        if not path.endswith(".html") or path.endswith("/index.html"):
            continue
        raw = git_show(repo, ref, path)
        txt = html.unescape(re.sub(r"<[^>]+>", " ", raw))
        txt = re.sub(r"\s+", " ", txt)
        for name in exports:
            if re.search(r"\b" + re.escape(name) + r"\s*\(", txt):
                covered.add(name)
    return covered


def check_worktree(repo: Path) -> tuple[list[str], set[str], set[str]]:
    ns = (repo / "NAMESPACE").read_text(encoding="utf-8")
    exports = exports_from_namespace_text(ns)
    aliases: set[str] = set()
    man_dir = repo / "man"
    for path in man_dir.glob("*.Rd"):
        aliases |= rd_aliases_from_text(path.read_text(encoding="utf-8"))
    rendered = rendered_coverage(repo / "docs" / "reference", exports)
    return exports, aliases, rendered


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--branch", default=None, help="Git ref to inspect (default: worktree)")
    parser.add_argument("--fail-if-missing", action="store_true")
    args = parser.parse_args()
    repo = args.repo.resolve()

    if args.branch:
        ns = git_show(repo, args.branch, "NAMESPACE")
        exports = exports_from_namespace_text(ns)
        aliases = collect_man_aliases(repo, args.branch)
        rendered = collect_rendered(repo, args.branch, exports)
        label = args.branch
    else:
        exports, aliases, rendered = check_worktree(repo)
        label = subprocess.check_output(
            ["git", "branch", "--show-current"], cwd=repo, text=True
        ).strip()

    missing_man = [e for e in exports if e not in aliases]
    missing_render = [e for e in exports if e not in rendered]

    print(f"branch/worktree: {label}")
    print(f"exports: {len(exports)}")
    print(f"man aliases covered: {len(exports) - len(missing_man)}/{len(exports)}")
    print(f"rendered reference covered: {len(rendered)}/{len(exports)}")
    if missing_man:
        print("\nMissing from man/*.Rd:")
        for name in missing_man:
            print(f"  - {name}")
    if missing_render:
        print("\nMissing from docs/reference:")
        for name in sorted(missing_render):
            print(f"  - {name}")

    if args.fail_if_missing and (missing_man or missing_render):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
