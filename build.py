"""
TODO:
- インデックスの生成
"""

import argparse
import os
import re
import subprocess
import time
import yaml

from pathlib import Path


def gather_markdown(paths: list[Path]):
    if len(paths) >= 1:
        return paths

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


def getMathjaxRelativePath(html: Path, online=False):
    if online:
        return ""
    mathjax_path = Path("./lib/MathJax/es5/tex-chtml-full.js")
    mathjax_rel_path = os.path.relpath(mathjax_path.resolve(), html.resolve())
    return "=" + "/".join(Path(mathjax_rel_path).parts[1:])


def pandoc_md_to_html5(md: Path, template_path: Path, rebuild=False, online=False):
    if md.suffix != ".md":
        return

    if str(md.name).lower() == "readme.md":
        return

    html = md.parent / Path(md.stem + ".html")
    if not is_source_modified(md, html) and not rebuild:
        return

    print("Processing " + str(md))

    result = subprocess.run(
        [
            "pandoc",
            "--lua-filter",
            "./pandocfilter.lua",
            "--standalone",
            "--toc",
            "--toc-depth=6",
            "--metadata",
            f"title={md.stem}",
            "--metadata",
            f"date={time.strftime('%F')}",
            "--metadata",
            "lang=ja",
            "--syntax-highlighting=some.theme",
            f"--template={str(template_path)}",
            "--from=markdown",
            "--to=html5",
            f"--mathjax{getMathjaxRelativePath(html, online)}",
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
    parser.add_argument(
        "-r", "--rebuild", action="store_true", help="Rebuild all files."
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Specify a single input markdown. The specified markdown is forced to be rebuild.",
        type=Path,
        nargs="*",
        default=[],
    )
    parser.add_argument(
        "--online",
        action="store_true",
        help="Change MathJax link to default CDN link provided by pandoc.",
    )
    args = parser.parse_args()

    if len(args.input) >= 1:
        args.rebuild = True

    dump_config_yml()

    mds = gather_markdown(args.input)
    for md in mds:
        pandoc_md_to_html5(md, "template.html", args.rebuild, args.online)
