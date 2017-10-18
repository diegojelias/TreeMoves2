import readTree

myTreeString = "(P:0.09,(Q:0.07,(X:0.02,((Y:0.03,Z:0.01):0.02,W:0.08):0.06):0.03):0.04);"

myTree = readTree.Tree(myTreeString)

myTree.print_names(myTree.root)

myTree.list_term_nodes(myTree.root)

myTree.tree_len(myTree.root)

myTreeDict = myTree.node_dict(myTree.root)

newTree = readTree.NNI(myTree)
newTree.newick(newTree.root)