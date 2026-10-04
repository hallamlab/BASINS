"""Render a code-based software workflow in the shared Hallam documentation style.

The JSON describes public software behavior, not manuscript results. Module rows
are conceptual groups; branching lanes represent evidence that converges.
"""
import html
import json


def render(root, stem="main-workflow", output="workflow-main"):
    data = json.loads((root / f'docs/diagrams/{stem}.json').read_text())
    width = 1720
    heights = [290 if 'lanes' in row else 190 for row in data['rows']]
    height = 310 + sum(heights)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
             f'<title id="title">{data["title"]} complete software workflow</title>',
             '<desc id="desc">Numbered conceptual modules with inputs, computational steps, data products and optional branches. See the workflow guide for source-code references and execution order.</desc>',
             '<rect width="100%" height="100%" fill="#ffffff"/>',
             '<defs><marker id="arrow" viewBox="0 0 10 10" refX="10" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0 0 L10 5 L0 10 Z" fill="#111111"/></marker></defs>',
             '<g font-family="Times New Roman, Times, serif" fill="#111111">']

    def text(x, y, label, size=22, anchor='middle'):
        for i, line in enumerate(label.split('\n')):
            parts.append(f'<text x="{x}" y="{y+i*(size+3)}" font-size="{size}" text-anchor="{anchor}">{html.escape(line)}</text>')

    def rect(x,y,w,h,fill,stroke,rx=0):
        parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="1.5"/>')

    def wire(points, arrow=True):
        path='M'+' L'.join(f'{x} {y}' for x,y in points)
        marker=' marker-end="url(#arrow)"' if arrow else ''
        parts.append(f'<path d="{path}" fill="none" stroke="#111111" stroke-width="1.5"{marker}/>')

    def symbol(x,y,kind):
        if kind=='compute':
            parts.append(f'<path d="M{x-13} {y} L{x} {y-13} L{x+13} {y} L{x} {y+13} Z" fill="#f5f5f5" stroke="#666666" stroke-width="3"/>')
        else:
            fill,stroke={'input':('#DAE8FC','#6C8EBF'),'output':('#D5E8D4','#82B366'),'data':('#F5F5F5','#666666')}[kind]
            parts.append(f'<circle cx="{x}" cy="{y}" r="12" fill="{fill}" stroke="{stroke}" stroke-width="3"/>')

    def node(x,y,n):
        symbol(x,y,n['kind'])
        text(x,y-25-(len(n['label'].split('\n'))-1)*25,n['label'])
        if n['tool']:text(x,y+37,n['tool'],18)

    # Title first, then inputs and the symbol key; module geometry stays fixed.
    rect(20,20,1680,78,'#CCCCCC','#666666',16)
    text(860,53,data['controller'],28)
    text(860,82,data['execution'],20)
    rect(20,120,975,130,'#DAE8FC','#6C8EBF',16)
    text(45,151,'Inputs',28,'start')
    for i,line in enumerate(data['inputs']):text(45,187+i*29,line,21,'start')
    rect(1030,120,670,110,'#CCCCCC','#666666',16)
    for x,label,kind in [(1090,'Module',None),(1220,'Compute','compute'),(1350,'Data','data'),(1480,'Input','input'),(1610,'Output','output')]:
        text(x,152,label,21)
        if kind:symbol(x,189,kind)
        else:rect(x-12,177,24,24,'#F5F5F5','#111111')
    top=290
    previous_y = None
    for index,(row,span) in enumerate(zip(data['rows'],heights),1):
        y=top+span/2-15
        if previous_y is not None:
            wire([(326, previous_y + 14), (326, y - 14)])
        previous_y = y
        # Module names are explicitly wrapped to fit the same label column.
        import textwrap
        title='\n'.join(textwrap.wrap(row['title'],22))
        text(160,y-(len(title.split('\n'))-1)*15,title,27)
        rect(313,y-13,26,26,'#F5F5F5','#111111')
        text(326,y+7,str(index),20)
        if 'lanes' in row:
            ys=[y-64,y+64]
            wire([(339,y),(350,y)],False)
            for lane,ly in zip(row['lanes'],ys):
                wire([(350,y),(350,ly),(448,ly)])
                for x,n in zip([460,820,1160],lane):node(x,ly,n)
                wire([(473,ly),(807,ly)])
                wire([(833,ly),(1147,ly)])
                wire([(1173,ly),(1370,ly),(1370,y)],False)
            wire([(1370,y),(1548,y)])
            node(1560,y,dict(label=row['result'],tool='',kind='output'))
        else:
            nodes=row['nodes'];xs=[460+i*1100/(len(nodes)-1) for i in range(len(nodes))]
            wire([(339,y),(448,y)])
            for x,n in zip(xs,nodes):node(x,y,n)
            for a,b in zip(xs,xs[1:]):wire([(a+13,y),(b-13,y)])
        if row.get('note'):text(1010,top+span-12,row['note'],18)
        top+=span
    parts.extend(['</g>','</svg>'])
    source=root/f'docs/assets/{output}.svg'
    source.write_text('\n'.join(parts)+'\n')
    return source,width,height
