#!/usr/bin/env python3
"""
Render a Markdown report to a self-styled HTML page for publishing as an Artifact
(no <html>/<head>/<body>: the publisher adds the skeleton). Light/dark theme-aware.

Usage: tools/md_to_artifact.py REPORT.md OUT.html --title "Page title"
"""
import argparse
import re
import subprocess

CSS = """
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:ital,wght@0,400;0,500;0,600;1,400&display=swap">
<style>
:root { --bg:#f6f7f9; --surface:#ffffff; --ink:#15171b; --ink2:#555a62; --ink3:#7c818a; --rule:#dde0e5;
        --accent:#2a78d6; --codebg:#eef0f3;
        --sans:"IBM Plex Sans","Helvetica Neue",Arial,sans-serif; --mono:"IBM Plex Mono",Menlo,Consolas,monospace; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { color-scheme:dark;
        --bg:#121417; --surface:#1a1d21; --ink:#eceef1; --ink2:#a9aeb6; --ink3:#858a93; --rule:#2d3137;
        --accent:#3987e5; --codebg:#22262b; } }
:root[data-theme="dark"] { color-scheme:dark; --bg:#121417; --surface:#1a1d21; --ink:#eceef1; --ink2:#a9aeb6;
        --ink3:#858a93; --rule:#2d3137; --accent:#3987e5; --codebg:#22262b; }
html { background: var(--bg); }
body { background: var(--bg); color: var(--ink); font-family: var(--sans); font-size: 16px; line-height: 1.6;
       margin: 0; padding-inline: 20px; padding-block: 36px 64px; }
main { max-width: 760px; margin: 0 auto; }
h1 { font-size: 30px; line-height: 1.2; font-weight: 600; margin: 0 0 12px; text-wrap: balance; letter-spacing: -.01em; }
h2 { font-size: 20px; font-weight: 600; margin: 36px 0 10px; text-wrap: balance; }
h3 { font-size: 16px; font-weight: 600; margin: 24px 0 8px; }
p, li { max-width: 68ch; }
p { margin: 0 0 12px; }
ul { padding-left: 22px; margin: 0 0 12px; }
li { margin: 4px 0; }
em { color: var(--ink2); }
a { color: var(--accent); }
code { font-family: var(--mono); font-size: .88em; background: var(--codebg); padding: 1px 5px; border-radius: 4px; }
hr { border: 0; border-top: 1px solid var(--rule); margin: 32px 0; }
.tw { overflow-x: auto; margin: 8px 0 16px; }
table { border-collapse: collapse; width: 100%; font-size: 14px; background: var(--surface); border: 1px solid var(--rule); border-radius: 8px; }
th, td { text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--rule); vertical-align: top; }
th { font-weight: 500; color: var(--ink2); font-size: 13px; }
tbody tr:last-child td { border-bottom: none; }
td:first-child { font-family: var(--mono); font-variant-numeric: tabular-nums; color: var(--ink2); }
strong { font-weight: 600; }
</style>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("md"); ap.add_argument("out"); ap.add_argument("--title", required=True)
    a = ap.parse_args()
    body = subprocess.run(["pandoc", a.md, "-f", "gfm", "-t", "html5"], capture_output=True, text=True, check=True).stdout
    body = re.sub(r"<table>", '<div class="tw"><table>', body)
    body = re.sub(r"</table>", "</table></div>", body)
    html = f"<title>{a.title}</title>\n{CSS}\n<main>\n{body}\n</main>\n"
    open(a.out, "w").write(html)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
