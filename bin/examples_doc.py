#!/usr/bin/env python3
"""
Prepares example notebooks for Sphinx documentation.

Copies notebooks from the examples/ directory into docs/examples/,
ordering them according to any ordering hints found in examples/README.md
(if present), and regenerates docs/examples.rst.
"""
import os
import re
import shutil
import json


def _extract_order(readme_path):
    """Return list of notebook filenames in the order they appear in README."""
    if not os.path.exists(readme_path):
        return []
    text = open(readme_path, "r", encoding="utf-8").read()
    return re.findall(r"MontuPython-[0-9A-Za-z_\-]+\.ipynb", text)


def _title_from_notebook(path):
    """Extract a display title from the first markdown heading in a notebook."""
    try:
        nb = json.load(open(path, "r", encoding="utf-8"))
        for cell in nb.get("cells", []):
            if cell.get("cell_type") != "markdown":
                continue
            for line in cell.get("source", []):
                text = line.strip()
                if text.startswith("#"):
                    return re.sub(r"^#+\s*", "", text).strip()
        return ""
    except Exception:
        return ""


def main():
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    src_dir = os.path.join(root_dir, "examples")
    docs_examples_dir = os.path.join(root_dir, "docs", "examples")
    rst_file = os.path.join(root_dir, "docs", "examples.rst")
    readme_path = os.path.join(src_dir, "README.md")

    # Recreate the docs/examples directory
    if os.path.exists(docs_examples_dir):
        shutil.rmtree(docs_examples_dir)
    os.makedirs(docs_examples_dir)

    # Build map of all available notebooks
    notebook_map = {
        name: os.path.join(src_dir, name)
        for name in os.listdir(src_dir)
        if name.endswith(".ipynb")
    }

    # Apply README ordering first, then append remaining notebooks sorted
    order = _extract_order(readme_path)
    ordered_files = []
    for name in order:
        path = notebook_map.pop(name, None)
        if path:
            ordered_files.append(path)
    for name in sorted(notebook_map.keys()):
        ordered_files.append(notebook_map[name])

    # Copy notebooks and collect names
    notebooks = []
    for src in ordered_files:
        basename = os.path.basename(src)
        dest = os.path.join(docs_examples_dir, basename)
        shutil.copy2(src, dest)
        notebooks.append(basename)
        print(f"Copied {basename} to docs/examples/")

    # Write examples.rst with explicit (and unique) titles so the sidebar
    # does not appear to have duplicated examples when notebooks share headings.
    with open(rst_file, "w", encoding="utf-8") as rst:
        rst.write("Examples\n")
        rst.write("========\n\n")
        rst.write(".. toctree::\n")
        rst.write("   :maxdepth: 2\n")
        rst.write("   :caption: Examples\n\n")
        seen_titles = set()
        for nb in notebooks:
            name = os.path.splitext(nb)[0]
            src = os.path.join(src_dir, nb)
            title = _title_from_notebook(src)
            if not title:
                title = name.replace("MontuPython-", "").replace("-", " ")

            # Ensure title uniqueness in the sidebar.
            if title in seen_titles:
                title = f"{title} ({name})"
            seen_titles.add(title)

            rst.write(f"   {title} <examples/{name}>\n")

    print("Done. examples.rst regenerated.")


if __name__ == "__main__":
    main()
