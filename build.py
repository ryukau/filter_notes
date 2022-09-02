"""
TODO:
- インデックスの生成
"""

import argparse
import json
import os
import re
import subprocess
import time
import yaml

from pathlib import Path

def gather_markdown():
    mds = []
    for path in Path(".").glob("*"):
        if path.is_dir() and not re.match(r"^[\._]", str(path)[0]):
            for md in path.glob("*.md"):
                mds.append(md)
    return mds

def is_source_modified(md, html):
    if not html.exists():
        return True
    return os.path.getmtime(md) > os.path.getmtime(html)

def pandoc_md_to_html5(md, template_path, rebuild=False):
    if md.suffix != ".md":
        return

    if str(md).lower() == "readme.md":
        return

    html = md.parent / Path(md.stem + ".html")
    if not is_source_modified(md, html) and not rebuild:
        return

    print("Processing " + str(md))

    result = subprocess.run(
        [
            "pandoc",
            "--standalone",
            "--toc",
            "--toc-depth=4",
            "--metadata",
            f"title={md.stem}",
            "--metadata",
            f"date={time.strftime('%F')}",
            "--metadata",
            "lang=ja",
            "--highlight-style",
            "some.theme",
            f"--template={str(template_path)}",
            "--from=markdown",
            "--to=html5",
            "--mathjax",
            f"--output={md.with_suffix('.html')}",
            md,
        ],
        stdout=subprocess.PIPE,
        encoding="utf-8",
    )

def dump_config_yml():
    md_list = sorted([str(md.as_posix()) for md in Path(".").glob("**/*.md")])
    with open("_config.yml", "w") as outfile:
        yaml.dump({"exclude": md_list}, outfile, default_flow_style=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rebuild", action='store_true', help="rebuild all file")
    args = parser.parse_args()

    dump_config_yml()

    mds = gather_markdown()
    for md in mds:
        pandoc_md_to_html5(md, "template.html", args.rebuild)
