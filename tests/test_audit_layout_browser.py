"""Measure actual DT geometry in the production WGS UI with synthetic data only."""
from __future__ import annotations

import json
import os
from pathlib import Path
import time
from urllib.request import urlopen

from playwright.sync_api import sync_playwright

URL = os.environ.get("AUDIT_TEST_URL", "http://127.0.0.1:8765")
OUT = Path("test-artifacts/audit-layout")
OUT.mkdir(parents=True, exist_ok=True)
AUDIT = 'a[data-value="DNA calls and sample audit"]'
OTHER = 'a[data-value="Genotype comparison"]'

GEOMETRY = r"""() => {
  const rect = e => { const r = e.getBoundingClientRect();
    return {top: r.top, bottom: r.bottom, height: r.height, left: r.left, right: r.right}; };
  const visible = e => e.getClientRects().length && getComputedStyle(e).visibility !== 'hidden';
  const measure = id => {
    const e = document.getElementById(id);
    const wrapper = e.querySelector('.dataTables_wrapper, .dt-container');
    const controls = Array.from(e.querySelectorAll('.dataTables_info, .dataTables_paginate, .dt-info, .dt-paging')).filter(visible);
    const bodies = Array.from(e.querySelectorAll('.dataTables_scrollBody, .dt-scroll-body')).filter(visible);
    const next = e.nextElementSibling;
    return {output: rect(e), wrapper: rect(wrapper), controls: controls.map(rect),
      bodies: bodies.map(rect), next: next ? rect(next) : null,
      nextText: next ? next.textContent.trim() : null,
      flexShrink: getComputedStyle(e).flexShrink, flexGrow: getComputedStyle(e).flexGrow};
  };
  return {samples: measure('audit-samples'), genotypes: measure('audit-genotypes')};
}"""


def problems(geometry: dict) -> list[str]:
    issues = []
    for field, table in geometry.items():
        if table["output"]["height"] <= 0:
            issues.append(f"{field}: output has no height")
        if len(table["controls"]) < 2:
            issues.append(f"{field}: missing info/pagination controls")
        for child in [table["wrapper"], *table["controls"], *table["bodies"]]:
            if child["bottom"] > table["output"]["bottom"] + 2:
                issues.append(f"{field}: child overflows its output")
            if table["next"] and child["bottom"] > table["next"]["top"] + 2:
                issues.append(f"{field}: child overlaps the following heading/note")
    if not geometry["samples"]["nextText"].startswith("All DNA calls"):
        issues.append("Missing second audit section heading")
    if geometry["samples"]["output"]["bottom"] > geometry["genotypes"]["output"]["top"] + 2:
        issues.append("Audit table output boxes overlap")
    return issues


def ready(page) -> None:
    page.locator(AUDIT).click()
    page.wait_for_function("""() => window.jQuery && jQuery.fn && jQuery.fn.dataTable &&
      jQuery.fn.dataTable.isDataTable(document.querySelector('#audit-samples .dataTables_scrollBody table')) &&
      jQuery.fn.dataTable.isDataTable(document.querySelector('#audit-genotypes .dataTables_scrollBody table')) &&
      document.querySelectorAll('#audit-samples .dataTables_scrollBody tbody tr').length === 10 &&
      document.querySelectorAll('#audit-genotypes .dataTables_scrollBody tbody tr').length === 10""")
    page.wait_for_timeout(500)


def check(page, label: str, results: list) -> None:
    # Wait for DT's asynchronous server-side draws and bslib resize handlers.
    deadline = time.monotonic() + 8
    while True:
        geometry = page.evaluate(GEOMETRY)
        issues = problems(geometry)
        if not issues or time.monotonic() > deadline:
            break
        page.wait_for_timeout(150)
    results.append({"case": label, "geometry": geometry, "issues": issues})
    (OUT / "geometry.json").write_text(json.dumps(results, indent=2))
    if issues:
        page.screenshot(path=str(OUT / f"FAILED-{label}.png"), full_page=True)
        raise AssertionError(f"{label}: {issues}")
    print(f"PASS {label}", flush=True)


def set_length(page, size: int) -> None:
    for field in ("samples", "genotypes"):
        page.locator(f'#audit-{field} .dataTables_length select').select_option(str(size))
    page.wait_for_function("""count => ['samples','genotypes'].every(field =>
      document.querySelectorAll('#audit-' + field + ' .dataTables_scrollBody tbody tr').length === count)""",
      arg=min(size, 75))
    page.wait_for_timeout(250)


def main() -> None:
    for _ in range(90):
        try:
            with urlopen(URL, timeout=2) as response:
                if response.status == 200:
                    break
        except OSError:
            time.sleep(1)
    else:
        raise RuntimeError("Synthetic Shiny fixture did not start")
    results: list = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        page.set_default_timeout(20000)
        errors = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        # Control: the original page must exhibit the reported overflow. The
        # baseline omits only presentation CSS, not tables or production UI.
        page.goto(URL + "?baseline=1")
        ready(page)
        baseline = page.evaluate(GEOMETRY)
        baseline_issues = problems(baseline)
        print("BASELINE", json.dumps({"issues": baseline_issues, "geometry": baseline}), flush=True)
        (OUT / "baseline.json").write_text(json.dumps(baseline, indent=2))
        page.screenshot(path=str(OUT / "before.png"), full_page=True)
        assert baseline_issues, "Baseline did not reproduce overflow; review fixture rather than claiming a regression test"

        page.goto(URL)
        ready(page)
        for size in (10, 25, 50, 100):
            set_length(page, size)
            check(page, f"desktop-{size}-rows", results)
        set_length(page, 10)
        horizontal = page.evaluate("""() => ['samples','genotypes'].map(field => {
          const body=document.querySelector('#audit-'+field+' .dataTables_scrollBody');
          body.scrollLeft=body.scrollWidth;
          return {offset:body.scrollLeft, overflow:getComputedStyle(body).overflowX};
        })""")
        assert all(x["offset"] > 0 and x["overflow"] in ("auto", "scroll") for x in horizontal)
        check(page, "horizontal-scroll-keeps-all-columns", results)
        page.evaluate("""() => document.querySelectorAll('.dataTables_scrollBody').forEach(e=>e.scrollLeft=0)""")
        page.locator('#audit-samples .dataTables_paginate .next').click()
        page.wait_for_function("""() => document.querySelector('#audit-samples .dataTables_info').textContent.includes('11')""")
        check(page, "next-page", results)
        # Search filters both to one row, then zero rows, without hiding the
        # following heading or the other table's controls.
        for term in ("synthetic-003", "NO-SUCH-SYNTHETIC-SAMPLE"):
            for field in ("samples", "genotypes"):
                page.locator(f'#audit-{field} .dataTables_filter input').fill(term)
            expected = 0 if term.startswith("NO-SUCH") else 1
            page.wait_for_function("""count => ['samples','genotypes'].every(field =>
              jQuery('#audit-' + field + ' .dataTables_scrollBody table').DataTable().page.info().recordsDisplay === count)""",
              arg=expected)
            page.wait_for_timeout(150)
            check(page, "search-" + term, results)
        for field in ("samples", "genotypes"):
            page.locator(f'#audit-{field} .dataTables_filter input').fill("")
        page.wait_for_function("""() => ['samples','genotypes'].every(field =>
          jQuery('#audit-' + field + ' .dataTables_scrollBody table').DataTable().page.info().recordsDisplay === 75)""")
        page.locator(OTHER).click()
        page.locator(AUDIT).click()
        page.wait_for_timeout(300)
        check(page, "return-from-hidden-tab", results)
        for width, height in ((1280, 800), (768, 900)):
            page.set_viewport_size({"width": width, "height": height})
            page.wait_for_timeout(500)
            check(page, f"viewport-{width}", results)
        page.set_viewport_size({"width": 1440, "height": 1000})
        page.evaluate("document.documentElement.style.zoom = '2'")
        page.wait_for_timeout(500)
        check(page, "css-zoom-200-percent", results)
        page.screenshot(path=str(OUT / "after-zoom.png"), full_page=True)
        page.evaluate("document.documentElement.style.zoom = '1'")
        page.wait_for_timeout(500)
        page.screenshot(path=str(OUT / "after.png"), full_page=True)
        assert not errors, f"Browser JavaScript errors: {errors}"
        browser.close()
    print(f"AUDIT_BROWSER_LAYOUT_OK ({len(results)} cases)", flush=True)


if __name__ == "__main__":
    main()
