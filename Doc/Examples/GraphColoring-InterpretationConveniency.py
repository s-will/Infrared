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
#
# As a toy example for the use of Infrared, we model (a small instance) of graph coloring. The graph coloring model
# defines variables for the graph nodes, whose possible values are the colors. Furthermore, we define inequality constraints for every edge, to ensure that connected nodes receive different colors in valid solutions.
#
# This model is then extended by an objective function that counts the different colors in cycles of size four.
# From the model, we create an Infrared solver to find optimal solutions. 
#
# For instructional purposes, we furthermore show the dependency graph of the model and the tree decomposition used by the solver.

# +
import infrared as ir

# specify the graph coloring instance
numcolors = 4 # colors 0..3
numnodes = 9 # nodes 1..9
edges = [(1,5),(2,3),(2,5),(2,7),(3,6),(5,6),(5,7),(5,8),(5,9),(6,8),(6,9),(7,8)]

# recall that we extend the standard graph coloring problem by counting different colors in cycles
cycles = [(2,3,5,6),(2,5,7,8),(5,6,7,8)]
# -

# Setup the model and construct a solver. Report the treewidth.

# +
model = ir.Model()

# Note: to keep this example simple, we add a dummy variable
# with index 0. In this way, we can use the 1-based node indices
# as variable indices (without the need for index transformations)
model.add_variables(1,(0,0))

# Add one variable per node
model.add_variables(numnodes,numcolors)

# define NotEquals constraints in Python
ir.def_constraint_class('NotEquals',
    lambda i,j: [i,j], # signature and dependencies
    lambda x, y: x!=y)  # constraint semantics

# add one inequality constraint per edge 
model.add_constraints(NotEquals(i,j) for i,j in edges)

# define cardinality for sets of 4 variables
ir.def_function_class('Card', # define "cardinality"
    lambda i, j, k, l: [i,j,k,l], # signature, dependencies
    lambda x,y,z,w: len({x,y,z,w}))  # function evaluation

# add the cardinality function for every cycle
model.add_functions([Card(i,j,k,l) for i,j,k,l in cycles], 'card')

model.set_feature_weight(-1,'card') # set up in order to *minimize*

solver = ir.Optimizer(model)
print(f'Tree width: {solver.treewidth()}')
# -

# ## Draw the dependency graph of the model
#
# Please ignore the dummy node 0 - we don't suppress it to keep things simple. Otherwise, the dependency graph is identical to the input graph as specified above.

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

# ## Plot tree decomposition
#
# This shows the tree decomposition of the dependency graph that is internally generated by the solver as basis for the efficient computation.

filename="treedecomp"
solver.plot_td(filename,'png')
Image(filename=filename+".png",width=400)


# ## Assignment interpretation
# ### Conversion of assignments to coloring dictionaries

def assignment_to_coloring(a):
        colors=["red","cyan","yellow","orange","magenta","green","blue","grey"]
        coloring = {i:colors[v] for i,v in enumerate(a.values())}
        del coloring[0]
        return coloring


# ### Syntactic sugar for more convenient use of assignments

# +
class ColoringAssignment(dict):
    def __init__(self, a):
        self.update(assignment_to_coloring(a))
        self._assignment = a

    def values(self):
        return self._assignment.values()
    
class ColoringOptimizer(ir.Optimizer):
    def __init__(self, model):
        super().__init__(model)
    
    def optimize(self):
        return ColoringAssignment(super().optimize())

class ColoringSampler(ir.Sampler):
    def __init__(self, model):
        super().__init__(model)
    
    def sample(self):
        return ColoringAssignment(super().sample())


# -

# ## Retreive colorings|

# +
model.set_feature_weight(-1,'card') # set negative weight to minimize card
solver = ColoringOptimizer(model)
best_coloring = solver.optimize()

#print('Optimal coloring (internal rep):', coloring.values())
print('Optimal coloring:', best_coloring)
print('Colors in cycles:', model.eval_feature(best_coloring,'card'))

# +
model.set_feature_weight(-1,'card') # set weight for number of colors in cycles
solver = ColoringSampler(model)
coloring = solver.sample()

print('Sampled coloring:', coloring)
print('Colors in cycles:', model.eval_feature(coloring,'card'))


# -

# ## Plot colorings

def plot_coloring(coloring):
    import graphviz

    G = graphviz.Graph('coloring', engine="fdp")

    for i in range(1,numnodes+1):
        G.node(f'{i}', style="filled", color=coloring[i])
    for i,j in edges:
        G.edge(f'{i}',f'{j}')

    return G


plot_coloring(best_coloring)

plot_coloring(coloring)


