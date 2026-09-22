"""Prepare a release: bump the version in Project.toml and close the NEWS section.

Mechanical half of docs/dev-guides/code-quality/changelogs_and_versions.md Part 3.
Registration stays a maintainer action: this only opens the PR.
"""

import argparse
import re
import sys

VERSION_RE = re.compile(r'^(version\s*=\s*")([^"]+)(")$', re.M)


def next_version(current, bump):
    major, minor, patch = (int(p) for p in current.split("."))
    # Julia's modified SemVer: below 1.0 the MINOR slot is the breaking slot,
    # so everything non-breaking lands in PATCH
    if major == 0:
        return f"0.{minor + 1}.0" if bump == "major" else f"0.{minor}.{patch + 1}"
    if bump == "major":
        return f"{major + 1}.0.0"
    if bump == "minor":
        return f"{major}.{minor + 1}.0"
    return f"{major}.{minor}.{patch + 1}"


def close_news_section(news, version):
    """Rename the open `main` section to `version` and open a fresh one."""
    lines = news.split("\n")

    def is_section_header(idx):
        # A header is a non-empty line underlined with dashes; the document
        # title is underlined with `=`, so it is not one
        return (
            idx + 1 < len(lines)
            and lines[idx].strip()
            and lines[idx + 1].strip()
            and set(lines[idx + 1].strip()) == {"-"}
        )

    for i, line in enumerate(lines):
        if line.strip() == "main" and is_section_header(i):
            break
    else:
        sys.exit("NEWS.md has no open 'main' section; see changelogs_and_versions.md 1.3")
    # The open section runs to the next header, i.e. the previous release
    end = next(
        (k for k in range(i + 2, len(lines)) if is_section_header(k)), len(lines)
    )
    released = "\n".join(lines[i + 2 : end]).strip()
    if not released:
        sys.exit("NEWS.md's 'main' section is empty: nothing to release")
    underline = lines[i + 1]
    lines[i : i + 2] = ["main", underline, "", version, "-" * len(version)]
    return "\n".join(lines), released


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bump", default="patch", choices=["patch", "minor", "major"])
    ap.add_argument("--version", default="")
    args = ap.parse_args()
    project = open("Project.toml").read()
    match = VERSION_RE.search(project)
    if not match:
        sys.exit("no version string in Project.toml")
    current = match.group(2)
    version = args.version.lstrip("v") or next_version(current, args.bump)
    news, released = close_news_section(open("NEWS.md").read(), f"v{version}")
    open("Project.toml", "w").write(VERSION_RE.sub(rf"\g<1>{version}\g<3>", project, count=1))
    open("NEWS.md", "w").write(news)
    print(f"current={current}")
    print(f"version={version}")
    print(f"released<<EOF\n{released}\nEOF")


if __name__ == "__main__":
    main()
