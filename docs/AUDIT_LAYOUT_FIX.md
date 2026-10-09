# Audit table layout regression

The DNA calls and sample audit page contains two independent DT tables. Each
must retain enough height for its rows, horizontal scrollbar, information and
pagination controls before the next heading begins.

## Change

`www/frontend.css` opts the namespaced `samples` and `genotypes` output containers
out of flex growing/shrinking and keeps their height intrinsic. A block formatting
context also contains DT's floated toolbar/footer. A separator marks the second
section. This is scoped to these two outputs, including new dataset namespaces;
it does not impose fixed row heights or change global card/plot sizing.

No table columns, values, filtering, pagination, exports, scientific calculations,
source access, job submission or scheduler settings are changed. Long file paths
still use DT's horizontal scrolling instead of being removed or truncated.

The browser regression fixture uses the production `variant_ui()` and extracts
the exact production audit-table renderer from `variant_server()`. Only the
result provider is replaced with 75 synthetic rows and synthetic paths. It does
not source `app.R`, connect a dataset, open a BAM/VCF or submit a job.

## Checks

The branch workflow starts this fixture on loopback and runs Chromium through
Playwright. It measures the output, table, footer and next heading rectangles,
first without the CSS fix (reproduction control), then with it. Cases include
10/25/50/100 rows, pagination, one/zero search matches, tab switching, laptop/narrow
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
