# Audit table layout regression

The DNA calls and sample audit page contains two independent DT tables. Each
must retain enough height for its rows, horizontal scrollbar, information and
pagination controls before the next heading begins.

## Change

The captured browser styles identified DT's rule
`.html-fill-container > .html-fill-item.datatables { flex-basis: 400px; }`.
The first CSS attempt disabled growing/shrinking but lost the basis to this
more-specific dependency rule, still leaving each output at 400px while its
rows and footer extended beyond it. The regression test caught this.

`www/frontend.css` explicitly overrides the complete flex shorthand for the
namespaced `samples` and `genotypes` outputs, including that basis, and keeps
their height intrinsic. A block formatting
context also contains DT's floated toolbar/footer. A separator marks the second
section. This is scoped to these two outputs, including new dataset namespaces;
it does not impose fixed row heights or change global card/plot sizing.

No table columns, values, filtering, pagination, exports, scientific calculations,
source access, job submission or scheduler settings are changed. Long file paths
remain on one line and use DT's horizontal scrolling instead of being removed
or truncated. This avoids very tall rows when a path contains many hyphens.

The browser regression fixture uses the production `variant_ui()` and extracts
the exact production audit-table renderer from `variant_server()`. Only the
result provider is replaced with 75 synthetic rows and synthetic paths. It does
not source `app.R`, connect a dataset, open a BAM/VCF or submit a job.

## Checks

The branch workflow starts this fixture on loopback and runs Chromium through
Playwright. It measures the output, table, footer and next heading rectangles,
first without the CSS fix (reproduction control), then with it. Cases include
10/25/50/100 rows, pagination, horizontal access to long columns, one/zero
server-confirmed search matches, tab switching, laptop/narrow
viewports and 2x CSS zoom. CSS zoom is a layout stress test, not a claim of testing
all native browser zoom implementations. Synthetic screenshots, geometry JSON
and the fixture log are saved as a short-lived workflow artifact, not committed
as participant data.

Run locally from the repository root with the R UI dependencies installed:

```sh
python -m pip install playwright==1.58.0
python -m playwright install chromium
Rscript tests/audit_layout_app.R
# In another terminal:
python tests/test_audit_layout_browser.py
```

The supported SCC session still needs a visual smoke check after updating and
restarting the app. Synthetic browser checks do not certify its installed package
versions or access to controlled data. Keep real screenshots and identifiers out
of the public repository.

References: bslib's "Filling layouts" article and DT's `dataTableOutput()` source
explain fill-item behavior and the DT-specific sizing contract.
