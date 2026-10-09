"""Capture synthetic computed styles and a script-free snapshot for layout debugging."""
import json
import time
from urllib.request import urlopen
from playwright.sync_api import sync_playwright
from test_audit_layout_browser import OUT, URL, ready, GEOMETRY

for _ in range(90):
    try:
        with urlopen(URL, timeout=2) as response:
            if response.status == 200:
                break
    except OSError:
        time.sleep(1)
else:
    raise RuntimeError("Synthetic Shiny fixture did not start")

with sync_playwright() as playwright:
    browser = playwright.chromium.launch()
    page = browser.new_page(viewport={"width": 1440, "height": 1000})
    page.goto(URL)
    ready(page)
    diagnostic = page.evaluate("""() => {
      const describe = e => {
        const s = getComputedStyle(e), r = e.getBoundingClientRect();
        const keys = ['display','position','height','minHeight','maxHeight','blockSize',
          'minBlockSize','maxBlockSize','flex','flexBasis','contain','overflowY','contentVisibility'];
        return {tag: e.tagName, id: e.id, className: e.className,
          inline: e.getAttribute('style'), rect: {top:r.top,bottom:r.bottom,height:r.height},
          computed: Object.fromEntries(keys.map(k=>[k,s[k]]))};
      };
      return ['audit-samples','audit-genotypes'].map(id=>{
        const e=document.getElementById(id);
        return {output:describe(e), parent:describe(e.parentElement),
          children:Array.from(e.children).map(describe),
          wrapper:describe(e.querySelector('.dataTables_wrapper,.dt-container'))};
      });
    }""")
    (OUT / "computed-styles.json").write_text(json.dumps(diagnostic, indent=2))
    snapshot = page.evaluate("""() => {
      let css = '';
      for (const sheet of document.styleSheets) {
        try { css += Array.from(sheet.cssRules).map(r=>r.cssText).join('\\n') + '\\n'; }
        catch (err) { css += '/* external stylesheet omitted */\\n'; }
      }
      const doc=document.documentElement.cloneNode(true);
      doc.querySelectorAll('script,link[rel=stylesheet]').forEach(e=>e.remove());
      const style=document.createElement('style'); style.textContent=css;
      doc.querySelector('head').appendChild(style);
      return '<!doctype html>\\n'+doc.outerHTML;
    }""")
    (OUT / "synthetic-layout-snapshot.html").write_text(snapshot)
    print("AUDIT_STYLE_DIAGNOSTICS", json.dumps(diagnostic), flush=True)
    browser.close()
