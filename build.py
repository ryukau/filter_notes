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

def get_match(md_info, md_name):
    for info in md_info:
        if info["name"] == md_name:
            return info
    return None

def read_build_info(rebuild=False):
    if not rebuild and Path("build_info").exists():
        with open("build_info", "r") as build_info:
            return json.load(build_info)
    return []

def gather_markdown():
    mds = []
    for path in Path(".").glob("*"):
        if path.is_dir() and not re.match(r"^[\._]", str(path)[0]):
            for md in path.glob("*.md"):
                mds.append(md)
    return mds

def get_last_modified(md):
    result = subprocess.run(
        ["git", "log", "-1", '--format="%ai"', "--", md],
        stdout=subprocess.PIPE,
        encoding="utf-8",
    )
    if len(result.stdout) > 0:
        return result.stdout[1:-1]  # remove double quotation (").

    # if there are no commits yet, get it from system.
    epoch = Path(md).stat().st_mtime
    return time.strftime("%Y-%m-%d %H:%M:%S %z", time.localtime(epoch))

def pandoc_md_to_html5(md, md_info, template_path):
    if md.suffix != ".md":
        return

    if str(md).lower() == "readme.md":
        return

    md_name = str(md)
    mtime = os.path.getmtime(md)
    info = get_match(md_info, md_name)
    if info is None:
        md_info.append({"name": md_name, "mtime": mtime})
    elif info["mtime"] == mtime:
        return
    else:
        info["mtime"] = mtime

    print("Processing " + md_name)

    result = subprocess.run(
        [
            "pandoc",
            "--standalone",
            "--toc",
            "--toc-depth=4",
            "--metadata",
            f"title={md.stem}",
            "--metadata",
            f"date={get_last_modified(md).split(' ')[0]}",
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
    md_list = sorted([str(md) for md in Path(".").glob("**/*.md")])
    with open("_config.yml", "w") as outfile:
        yaml.dump({"exclude": md_list}, outfile, default_flow_style=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rebuild", action='store_true', help="rebuild all file")
    args = parser.parse_args()

    dump_config_yml()

    md_info = read_build_info(args.rebuild)
    mds = gather_markdown()
    for md in mds:
        pandoc_md_to_html5(md, md_info, "template.html")

    with open("build_info", "w") as build_info:
        json.dump(md_info, build_info)
