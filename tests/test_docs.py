"""The documentation must not drift from the repository: every path the
README and CONTRIBUTING tables name exists (or is declared generated), and
every command AGENTS.md and CONTRIBUTING.md tell a contributor to run names
a script that exists."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENERATED = {"demo/"}            # documented as generated and not committed


def _table_paths(name):
    text = (ROOT / name).read_text(encoding="utf-8")
    paths = set()
    for line in text.splitlines():
        if line.startswith("| `"):
            first_cell = line.split("|")[1]
            paths |= set(re.findall(r"`([^`]+)`", first_cell))
    return paths


def test_every_path_in_the_file_tables_exists():
    for doc in ("README.md", "CONTRIBUTING.md"):
        for p in _table_paths(doc):
            if p in GENERATED:
                continue
            assert (ROOT / p.rstrip("/")).exists(), f"{doc} names {p}, which does not exist"


def test_every_documented_command_names_an_existing_script():
    for doc in ("AGENTS.md", "CONTRIBUTING.md"):
        text = (ROOT / doc).read_text(encoding="utf-8")
        for script in re.findall(r"python (?!-m)([\w./]+\.py)", text):
            assert (ROOT / script).exists(), f"{doc} runs {script}, which does not exist"
