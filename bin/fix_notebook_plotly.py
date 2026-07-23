#!/usr/bin/env python3
import sys
import nbformat
import plotly.io as pio
import plotly.graph_objects as go

def fix_plotly_notebook(path):
    print(f"Fixing notebook for docs (Plotly HTML injection): {path}")
    nb = nbformat.read(path, as_version=4)
    modified = False
    
    for cell in nb.cells:
        if cell.cell_type == 'code':
            for out in cell.outputs:
                if 'data' in out and 'application/vnd.plotly.v1+json' in out['data']:
                    plotly_data = out['data']['application/vnd.plotly.v1+json']
                    
                    # Upgrade scattermapbox to scattermap (MapLibre) to avoid blank maps in Sphinx
                    import json
                    json_str = json.dumps(plotly_data)
                    json_str = json_str.replace('"scattermapbox"', '"scattermap"')
                    json_str = json_str.replace('"mapbox"', '"map"')
                    json_str = json_str.replace('"mapbox_style"', '"map_style"')
                    plotly_data = json.loads(json_str)

                    fig = go.Figure(plotly_data)
                    # Force overwrite with 'cdn' version to ensure it works on ReadTheDocs without require.js dependencies
                    html = pio.to_html(fig, include_plotlyjs='cdn', full_html=False)
                    out['data']['text/html'] = html
                    modified = True
                        
    if modified:
        nbformat.write(nb, path)
        print(f"Saved modified notebook: {path}")

if __name__ == '__main__':
    for path in sys.argv[1:]:
        fix_plotly_notebook(path)
