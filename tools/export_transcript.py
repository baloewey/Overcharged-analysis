#!/usr/bin/env python3
"""
Export a Claude Code session transcript to a readable PDF.

Includes: every user message and assistant reply (timestamps in Chile time, UTC-3), shell
commands the user ran with `!` (with output), and one grey line per tool call (tool + its
description). Leaves out: full tool outputs, background-task notifications, internal reasoning.

Pipeline: session .jsonl -> Markdown (GitHub-flavoured) -> HTML (pandoc) -> PDF (headless Chrome).
Needs: pandoc, google-chrome (both on the office machine).

Usage:
    tools/export_transcript.py                      # latest session of this project -> PDF in project root
    tools/export_transcript.py -o notes.pdf         # choose output path
    tools/export_transcript.py --session FILE.jsonl # a specific session
    tools/export_transcript.py --md-only -o x.md    # just the Markdown
"""
import argparse
import datetime
import glob
import json
import os
import re
import subprocess
import sys
import tempfile

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SESSIONS = os.path.expanduser("~/.claude/projects/-home-baloewe-Physics-Overcharged")
TZ = datetime.timezone(datetime.timedelta(hours=-3))

CSS = """
@page { size: A4; margin: 16mm 15mm; }
body { font-family: "DejaVu Sans", "Liberation Sans", sans-serif; font-size: 10pt; line-height: 1.45;
       color: #0b0b0b; max-width: none; }
h1.title { font-size: 18pt; margin-bottom: 0; }
p.subtitle { font-size: 11pt; color: #52514e; margin-top: 2px; }
p.date { color: #52514e; }
h2 { font-size: 11pt; margin: 18px 0 6px; padding: 4px 8px; border-radius: 4px; page-break-after: avoid; }
h2.user { background: #e8f1fc; color: #184f95; }
h2.claude { background: #f1f0ec; color: #3d3c39; }
h3 { font-size: 10pt; color: #52514e; }
.tool { font-family: "DejaVu Sans Mono", monospace; font-size: 7.5pt; color: #8a8983; margin: 0; line-height: 1.3; }
pre, code { font-family: "DejaVu Sans Mono", monospace; font-size: 8pt; }
pre { background: #f6f6f3; padding: 6px 8px; border-radius: 4px; white-space: pre-wrap; word-break: break-word; }
table { border-collapse: collapse; font-size: 8.5pt; margin: 6px 0; }
th, td { border: 1px solid #d9d8d2; padding: 2px 6px; }
th { background: #f1f0ec; }
blockquote { color: #52514e; border-left: 3px solid #d9d8d2; margin-left: 0; padding-left: 10px; }
"""


def ts(s):
    try:
        return datetime.datetime.fromisoformat(s.replace("Z", "+00:00")).astimezone(TZ)
    except Exception:
        return None


def tool_line(b):
    name = b.get("name", "tool")
    inp = b.get("input", {}) or {}
    if name in ("Write", "Edit", "Read") and inp.get("file_path"):
        desc = inp["file_path"].replace(os.path.expanduser("~") + "/", "~/")
    else:
        desc = (inp.get("description") or str(inp.get("prompt", ""))[:80] or inp.get("file_path")
                or inp.get("url") or inp.get("query") or inp.get("skill") or "")
    desc = re.sub(r"\s+", " ", str(desc)).strip().replace("<", "&lt;")
    return f'<div class="tool">▸ {name}: {desc}</div>'


def to_markdown(src):
    lines, pending, last_role = [], [], None
    t_first = t_last = None

    def flush():
        nonlocal pending
        if pending:
            lines.extend(pending + [""])
            pending = []

    for raw in open(src):
        d = json.loads(raw)
        typ = d.get("type")
        if typ not in ("user", "assistant") or d.get("isMeta"):
            continue
        c = (d.get("message") or {}).get("content")
        t = ts(d.get("timestamp", ""))
        if t:
            t_first = t_first or t
            t_last = t
        when = t.strftime("%Y-%m-%d %H:%M") if t else ""
        if typ == "user":
            texts = [c] if isinstance(c, str) else \
                [b.get("text", "") for b in (c or []) if isinstance(b, dict) and b.get("type") == "text"]
            for txt in texts:
                if "<task-notification>" in txt:
                    continue
                txt = re.sub(r"<system-reminder>.*?</system-reminder>", "", txt, flags=re.S).strip()
                if not txt:
                    continue
                bi = re.search(r"<bash-input>(.*?)</bash-input>", txt, re.S)
                if bi:
                    so = re.search(r"<bash-stdout>(.*?)</bash-stdout>", txt, re.S)
                    se = re.search(r"<bash-stderr>(.*?)</bash-stderr>", txt, re.S)
                    outp = ((so.group(1) if so else "") + (se.group(1) if se else "")).strip()
                    flush()
                    lines += [f"### User ran a shell command · {when}", "", "```",
                              "$ " + bi.group(1).strip(), outp[:3000], "```", ""]
                    continue
                if "<bash-stdout>" in txt or "<bash-stderr>" in txt:
                    continue
                flush()
                lines += [f"## User · {when}", "", txt, ""]
                last_role = "user"
        elif isinstance(c, list):
            for b in c:
                if b.get("type") == "tool_use":
                    pending.append(tool_line(b))
                elif b.get("type") == "text" and b.get("text", "").strip():
                    flush()
                    if last_role != "assistant":
                        lines += [f"## Claude · {when}", ""]
                    lines += [b["text"].strip(), ""]
                    last_role = "assistant"
    flush()
    span = ""
    if t_first and t_last:
        span = t_first.strftime("%Y-%m-%d") + (
            "" if t_first.date() == t_last.date() else " to " + t_last.strftime("%Y-%m-%d"))
    return "\n".join(lines) + "\n", span, t_first, t_last


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--session", help="session .jsonl (default: most recently modified in this project)")
    ap.add_argument("-o", "--out", help="output file (default: PROJECT/Overcharged_session_<dates>.pdf)")
    ap.add_argument("--md-only", action="store_true", help="write Markdown instead of PDF")
    a = ap.parse_args()

    src = a.session or max(glob.glob(os.path.join(SESSIONS, "*.jsonl")), key=os.path.getmtime)
    md, span, t0, t1 = to_markdown(src)
    if a.out:
        out = a.out
    else:
        tag = t0.strftime("%Y-%m-%d") + ("" if t0.date() == t1.date() else "_" + t1.strftime("%d"))
        out = os.path.join(PROJECT, f"Overcharged_session_{tag}." + ("md" if a.md_only else "pdf"))
    if a.md_only:
        open(out, "w").write(md)
        print("wrote", out)
        return

    with tempfile.TemporaryDirectory() as tmp:
        mdf, css, html = (os.path.join(tmp, n) for n in ("t.md", "style.css", "t.html"))
        open(mdf, "w").write(md)
        open(css, "w").write(CSS)
        subprocess.run(["pandoc", mdf, "-f", "gfm+tex_math_dollars", "-t", "html5", "--standalone",
                        "--embed-resources", "--css", css, "--mathjax",
                        "--metadata", "title=Overcharged — working session with Claude Code",
                        "--metadata", "subtitle=Active nematic MPCD: defects in cavities with custom anchoring",
                        "--metadata", f"date={span}", "-o", html],
                       check=True, stderr=subprocess.DEVNULL)
        s = open(html).read()
        s = re.sub(r'<h2 id="([^"]*)">User', r'<h2 id="\1" class="user">User', s)
        s = re.sub(r'<h2 id="([^"]*)">Claude', r'<h2 id="\1" class="claude">Claude', s)
        open(html, "w").write(s)
        r = subprocess.run(["google-chrome", "--headless=new", "--disable-gpu", "--no-sandbox",
                            "--no-pdf-header-footer", f"--print-to-pdf={os.path.abspath(out)}",
                            "file://" + html], capture_output=True, text=True)
        if not os.path.exists(out):
            sys.exit("PDF not written:\n" + r.stderr[-2000:])
    print("wrote", out)


if __name__ == "__main__":
    main()
