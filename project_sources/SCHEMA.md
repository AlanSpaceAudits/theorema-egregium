# project_sources layout

- `pdf/`      source files (PDF, scans, OCR text)
- `web/`      saved copies of web pages (PDF print or HTML)
- `shots/`    page images of the exact passages cited in the deck
- `railway/`  everything for the Trans-Australian Railway distance, same layout inside
- `citations.json` and `railway/citations.json`: one entry per source

## citations.json entry

```
"<key>": {
  "apa": "APA 7 reference",
  "url": "where it was fetched from, or null",
  "file": "path relative to project_sources, or null",
  "license": "public domain | copyrighted, excerpt only | ...",
  "shots": [
    { "img": "shots/<key>_<what>.png", "page": "p. 20", "caption": "one line: what this image shows",
      "quote": "verbatim text of the passage the deck relies on" }
  ]
}
```

Rules: quotes are verbatim from the source. No invented page numbers, figures or sources. If a source cannot be
obtained, the entry still exists with `"file": null`, `"shots": []` and a `"missing": "what was tried"` field.
Page images: PNG, about 1400 px wide, the cited passage marked with a coloured box where practical.

## translations/

Every cited page in a language other than English has, in `translations/`:

- `<key>_pNNN_orig.png`  the full original page
- `<key>_pNNN_orig.txt`  a transcription of the original text
- `<key>_pNNN_en.tex`    the page recreated in English, in LaTeX, laid out to match the original
- `<key>_pNNN_en.pdf` and `_en.png`  the compiled English page

`translations_<lang>.json` lists them. The build shows the original and the English page side by side in the
sources reader. A `checked` field records the independent check of the translation.
