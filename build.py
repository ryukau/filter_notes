"""
TODO:
- タイトルの書き換え
- インデックスの生成
"""

import argparse
import json
import os
import re
import subprocess
import yaml

from bs4 import BeautifulSoup
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
        if path.is_dir() and not re.match("^[\._]", str(path)[0]):
            for md in path.glob("*.md"):
                mds.append(md)
    return mds

def pandoc_md_to_html5(md, md_info):
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
        ["pandoc", "-f", "markdown", "-t", "html5", "--katex", md],
        stdout=subprocess.PIPE,
        encoding="utf-8")

    html = template_html.replace("<!-- replace me -->", result.stdout)

    with open(md.with_suffix(".html"), "w") as html_file:
        print(html, file=html_file)

def dump_config_yml():
    md_list = sorted([str(md) for md in Path(".").glob("**/*.md")])
    print()
    with open("_config.yml", "w") as outfile:
        yaml.dump({"exclude": md_list}, outfile, default_flow_style=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r", "--rebuild", action='store_true', help="rebuild add file")
    args = parser.parse_args()

    dump_config_yml()

    with open("template.html", "r") as temp:
        template_html = temp.read()
    md_info = read_build_info(args.rebuild)
    mds = gather_markdown()
    for md in mds:
        pandoc_md_to_html5(md, md_info)

    with open("build_info", "w") as build_info:
        json.dump(md_info, build_info)
