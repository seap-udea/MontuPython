"""Build compact Plotly HTML for Qt WebEngine (no inline JS bundle)."""

from __future__ import annotations

from pathlib import Path

import plotly.graph_objects as go


def plotly_js_path() -> Path:
    import plotly

    return Path(plotly.__file__).parent / "package_data" / "plotly.min.js"


_PLOTLY_RESIZE_SCRIPT = """
<script>
function montuResizePlotly() {
  if (!window.Plotly) return;
  var minH = window.MONTU_PLOTLY_MIN_HEIGHT || 320;
  var viewH = window.innerHeight || document.documentElement.clientHeight || minH;
  var fill = window.MONTU_PLOTLY_FILL_VIEWPORT !== false;
  var h = fill ? Math.max(minH, viewH) : minH;
  document.body.style.minHeight = h + 'px';
  document.querySelectorAll('.plotly-graph-div').forEach(function(el) {
    el.style.height = h + 'px';
    el.style.minHeight = h + 'px';
    Plotly.Plots.resize(el);
  });
}
window.addEventListener('load', function() {
  montuResizePlotly();
  setTimeout(montuResizePlotly, 100);
  setTimeout(montuResizePlotly, 400);
});
window.addEventListener('resize', montuResizePlotly);
if (typeof ResizeObserver !== 'undefined') {
  new ResizeObserver(montuResizePlotly).observe(document.body);
}
</script>
"""


def figure_to_html(
    fig: go.Figure,
    *,
    min_height: int = 400,
    fill_viewport: bool = True,
    page_background: str = "#ffffff",
) -> str:
    """Return a small standalone HTML page for a Plotly figure.

    When *fill_viewport* is True (default), the plot height is
    ``max(min_height, viewport height)`` so tall panels are filled and
    short panels still scroll if *min_height* exceeds the view.
    """
    div = fig.to_html(include_plotlyjs=False, full_html=False, config={
        "responsive": True,
        "displayModeBar": True,
    })
    js_url = plotly_js_path().as_uri()
    h = max(320, int(min_height))
    fill_flag = "true" if fill_viewport else "false"
    return f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8"/>
  <script src="{js_url}"></script>
  <script>
    window.MONTU_PLOTLY_MIN_HEIGHT = {h};
    window.MONTU_PLOTLY_FILL_VIEWPORT = {fill_flag};
  </script>
  <style>
    html, body {{
      margin: 0;
      padding: 0;
      width: 100%;
      min-height: 100%;
      overflow: auto;
      background: {page_background};
    }}
    .plotly-graph-div {{
      width: 100% !important;
    }}
    /* Keep modebar off the title — park it above the rangeslider */
    .js-plotly-plot .plotly .modebar {{
      top: auto !important;
      left: auto !important;
      right: 12px !important;
      bottom: 52px !important;
    }}
    .js-plotly-plot .plotly .modebar-group {{
      background: rgba(255, 255, 255, 0.85) !important;
      border-radius: 4px;
      box-shadow: 0 1px 4px rgba(0, 0, 0, 0.12);
    }}
  </style>
</head>
<body>
{div}
{_PLOTLY_RESIZE_SCRIPT}
</body>
</html>"""
