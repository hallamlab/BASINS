"""Draw the public workflow SVG from its versioned JSON description."""
import html,json

def render(root):
    data=json.loads((root/'docs/diagrams/main-workflow.json').read_text())
    name,rows=data['title'],data['rows']
    width,height=1400,110+len(rows)*172
    colors=['#d9e5ef','#eadff0','#dceadb','#f5e5cb','#e3e3ef','#d6e9e7','#eedfdc','#e1e1e1']
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">',f'<title>{name} complete workflow</title>','<rect width="100%" height="100%" fill="white"/>','<defs><marker id="arrow" viewBox="0 0 10 10" refX="10" refY="5" markerWidth="8" markerHeight="8" orient="auto"><path d="M0 0 L10 5 L0 10 Z" fill="#666"/></marker></defs>',f'<text x="700" y="46" text-anchor="middle" font-family="Georgia,serif" font-size="32" font-weight="bold">{name} workflow</text>']
    for i,(a,b,c) in enumerate(rows):
        y=80+i*172
        parts += [f'<rect x="22" y="{y}" width="1356" height="142" rx="12" fill="{colors[i%len(colors)]}" stroke="#777" stroke-width="2"/>',f'<text x="48" y="{y+35}" font-family="Georgia,serif" font-size="25" font-weight="bold">{i+1}. {html.escape(a)}</text>',f'<text x="48" y="{y+77}" font-family="Georgia,serif" font-size="25">{html.escape(b)}</text>',f'<text x="48" y="{y+115}" font-family="Georgia,serif" font-size="22">{html.escape(c)}</text>']
        if i<len(rows)-1: parts += [f'<path d="M700 {y+142} V{y+172}" stroke="#666" stroke-width="2.5" fill="none" marker-end="url(#arrow)"/>']
    parts+=['</svg>']
    source=root/'docs/assets/workflow-main.svg'
    source.write_text('\n'.join(parts)+'\n')
    return source,width,height
