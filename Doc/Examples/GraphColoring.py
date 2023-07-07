# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Graph coloring

# +
import infrared as ir

numcolors = 4
#nodes = list(range(1,10))
#edges = [(1,5),(2,3),(2,5),(2,7),(3,6),(5,6),(5,7),(5,9),(6,8),(6,9),(7,8)]

nodes = list(range(1,16))
edges = [(1,5),(2,3),(2,5),(2,7),(3,6),(5,6),(5,7),(5,9),(6,8),(6,9),(7,8)]
for i in range(10,16):
    for j in range(i+1,16):
        edges.append((i,j))

cycles = [(2,3,5,6),(5,6,7,8),(2,5,7,8)]

nidx={x:i for i,x in enumerate(nodes)}
numnodes = len(nodes)

model = ir.Model(numnodes,numcolors)

ir.def_constraint_class('DiffColor',
    lambda i,j: [i,j],
    lambda x,y: x != y)

ir.def_function_class('UnusedNum',
    lambda i,j,k,l: [i,j,k,l],
    lambda x,y,z,w: len(set(range(4)-set([x,y,z,w]))))

model.add_constraints(DiffColor(nidx[i],nidx[j]) for i,j in edges)

model.add_functions(UnusedNum(*[nidx[i] for i in x]) for x in cycles)

solver = ir.Optimizer(model)
print(solver.treewidth())

# +
from IPython.display import Image
import re

# Plot dependency graph
filename = 'dependency_graph.dot'
model.write_graph(filename, True)

ir.dotfile_to_png(filename)
filename = re.sub(r"dot$","png",filename)

Image(filename=filename,width=300)
# -

# Plot tree decomposition
filename="treedecomp"
solver.plot_td(filename,'png')
Image(filename=filename+".png",width=400)


