"""index.html from README.md (run after editing the README): python3 make_page.py"""
import re
from pathlib import Path

import markdown

md = Path(__file__).with_name("README.md").read_text()
body = markdown.markdown(md, extensions=["tables", "fenced_code"])
body = re.sub(r'(?<!href=")(https?://[^\s<"]+)(?=[\s<])', r'<a href="\1">\1</a>', body)
for f in ("iq_to_trace.py", "fit_pass.py", "make_test_iq.py"):
    body = body.replace(f"<code>{f}</code>", f'<a href="{f}"><code>{f}</code></a>', 1)
css = """:root{--bg:#f7f8f8;--fg:#1b2422;--mut:#55625f;--rule:#d5dbd9;--acc:#2a5f9e;--code:#eceff0}
@media (prefers-color-scheme:dark){:root{--bg:#121716;--fg:#e3eae8;--mut:#a2b0ac;--rule:#2f3b39;--acc:#82b0e8;--code:#1f2928;color-scheme:dark}}
body{margin:0;background:var(--bg);color:var(--fg);font:17px/1.6 "Segoe UI",system-ui,sans-serif}
main{max-width:880px;margin:0 auto;padding:36px max(16px,env(safe-area-inset-left)) 64px}
h1{font-size:2.2rem;margin:0 0 .4em}h2{border-top:2px solid var(--fg);padding-top:10px;margin-top:2em}
a{color:var(--acc)}code{background:var(--code);padding:1px 5px;border-radius:3px;font-size:.88em}
pre{background:var(--code);padding:12px 14px;border-radius:6px;overflow-x:auto}pre code{background:none;padding:0}
table{border-collapse:collapse;width:100%;font-size:.94rem;display:block;overflow-x:auto}th,td{border-bottom:1px solid var(--rule);padding:6px 10px;text-align:left;vertical-align:top}
img{max-width:100%;height:auto;border:1px solid var(--rule);border-radius:4px;display:block;margin:10px 0}"""
page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Doppler Kit</title><style>{css}</style></head>
<body><main>{body}</main></body></html>
"""
Path(__file__).with_name("index.html").write_text(page)
print("index.html", len(page))
