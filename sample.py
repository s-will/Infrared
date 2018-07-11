import sample as sm


## Construct constraint network

var = { name:index for index,name in enumerate(["x","y","z"]) }

## for convenience (demonstrating for a simple case how more high
## level interfaces can be written on the Python side)
def binary_fun(x,y,f):
    return Function(lambda a: f(a[var[x]], a[var[x]]), [var[x],var[y]])

def binary_con(x,y,c):
    return Constraint(lambda a: c(a[var[x]], a[var[x]]), [var[x],var[y]])


cn = sm.ConstraintNetwork( [ range(0,4) ] * len(var) ) ## list of variable domains; list length = number of variables


vars_f1 = [var["x"],var["y"]]
f1 = Function( lambda a: a[var["x"]] + a[var["y"]], vars_f1 )


# next function def with little more sugare
f2 = binary_fun( "x", "z", lambda x,z: x+z )

vars_c1 = [var["x"],var["y"]]
c1 = Constraint( lambda a: a[var["x"]] < a[var["y"]], vars_c1 )

cn.add(f1)
cn.add(f2)
cn.add(c1)

## Construct tree decomposition
td = TreeDecomposition(cn)

## since the cn knows about functions and constraints associated to vars,
## we can support automatic assignment of functions and constraints to cluster,
## however note that there may be choices; we can choose to assign greedily
## in the order of cluster additions

cl1 = Cluster( [var["x"],var["y"]] )
td.add( cl1 )

cl2 = Cluster( [var["y"],var["z"]] )
td.addChildren( cl1, [ cl2 ] )


## Generate samples

## td.prepare_sampling()

for i in range(0:1000):
    sample = td.sample()


############################################################
## use of predefined functions and constraints

## complementary(x,y) --- beware of inefficiency due to dictionary indirection
complconstr = binary_con("x","y",complementary)
## better: (where ComplementConstraint is defined on C++ side to construct the appropriate Constraint object)
complconstr = ComplementConstraint(var["x"],var["y"])

stackEfunction = StackEnergy(i,j,k,l)
