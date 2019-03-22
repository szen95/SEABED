import sys
import graphviz as gv
from colour import Color
import modify_tree_image as mti

# the styles applied to the whole graphviz graph
STYLES = {
    'graph': {
        'splines': 'ortho',
        'overlap': 'scale',
        'nodesep' : '1',
        'rank': 'same',
    },
    'nodes': {
        'fontsize': '11',
        'margin': '0',
        'fontname': 'calibri',
        'height': '0.5',
        'width': '1',
    },
    'edges': {
        'rank': 'same',
    }
}

# applies the given style rules to a graphviz tree
def apply_styles(graph, styles):
    graph.graph_attr.update(
        ('graph' in styles and styles['graph']) or {}
    )
    graph.node_attr.update(
        ('nodes' in styles and styles['nodes']) or {}
    )
    graph.edge_attr.update(
        ('edges' in styles and styles['edges']) or {}
    )
    return graph

# a Node in the tree
class Node:
    def __init__(self, pos, value=None, outcome=None, parent=None):
        self.value = value
        self.pos = pos
        self.outcome = outcome
        self.parent = parent
        self.homo = False
        self.left = None
        self.right = None

    # checks if the given node is a leaf
    def is_leaf(self):
        return self.left  == None and self.right == None

    # gets the graphviz edge style
    def incoming_edge_style(self):
        return 'dotted' if self.homo else ''

    # sets the graphviz attributes for the node
    def set_visual_attributes(self, colors, xlab):
        # color
        self.color = 'gray'
        self.outcome = min(self.outcome, 99)
        if self.parent == None:
            self.color = 'red' if self.outcome == None else colors[self.outcome]
        if self.is_leaf():
            self.color = "gray"
            # self.color = "#4c4cff"
            print 'This is a leaf.'
            # perc = '.3' if len(xlab) == 1 else '.35'
            # self.color = "black;{}:cyan".format(perc) if self.outcome == None \
            #     else "black;{}:{}".format(perc, colors[self.outcome])
        # outline color
        self.outline = 'black'
        # shape
        self.shape = 'box' if self.is_leaf() else 'box'
        # label
        v = str(self.value if int(self.value) != self.value else int(self.value))
        if self.is_leaf():
            # self.label = '<<b><font color="white">  {}    </font>{}     </b>>'.format(xlab, v)
            self.label = v
        else:
            self.label = v

# gets the direction from the current node towards the given position (0 is left and 1 is right)
def direction(cur, pos):
    return int(pos.replace(cur.pos, '')[0])

# returns whether or not a node with the given position is in the tree
def in_tree(cur, pos):
    return closest(cur, pos).pos == pos

# finds the node in the tree that is closest to the given position
def closest(cur, pos):
    if cur.pos == pos:
        return cur
    nxt = direction(cur, pos)
    if nxt == 0:
        if cur.left == None:
            return cur
        else:
            return closest(cur.left, pos)
    else:
        if cur.right == None:
            return cur
        else:
            return closest(cur.right, pos)

# prints the tree to the screen
def print_tree(node, indent=0):
    print ' '*indent, node.pos, '|', node.value, \
        '| o='+str(node.outcome) if node.outcome != None else '', \
        '| ' + str(node.homo) if node.is_leaf() else ''
    if node.left != None:
        print_tree(node.left, indent + 1)
    if node.right != None:
        print_tree(node.right, indent + 1)

# sets the values of the nodes based on their children
def set_values(node, leaves_map, outcome_map, homo_map):
    if node.left == None and node.right == None:
        node.value = leaves_map[node.pos]
        if outcome_map != None:
            node.outcome = outcome_map[node.pos]
        node.homo = homo_map[node.pos]
    elif node.left == None:
        set_values(node.right, leaves_map, outcome_map, homo_map)
        node.value = node.right.value
    elif node.right == None:
        set_values(node.left, leaves_map, outcome_map, homo_map)
        node.value = node.left.value
    else:
        set_values(node.left, leaves_map, outcome_map, homo_map)
        set_values(node.right, leaves_map, outcome_map, homo_map)
        node.value = node.left.value + node.right.value

# handles the creation of the visual tree in graphviz
def draw_tree(graph, node, colors, leafNames):
    global cur_node_num
    global node_names
    xlab = str(leafNames.index(node.pos) + 1) if node.is_leaf() else ''
    node.set_visual_attributes(colors, xlab)
    node_names[node.pos] = str(cur_node_num)
    graph.node(str(cur_node_num), label=node.label, shape=node.shape,
        style='filled', fillcolor=node.color)
    cur_node_num += 1
    if node.left != None:
        left = draw_tree(graph, node.left, colors, leafNames)
        graph.edge(node_names[node.pos], node_names[left], style=node.left.incoming_edge_style())
    if node.right != None:
        right = draw_tree(graph, node.right, colors, leafNames)
        graph.edge(node_names[node.pos], node_names[right], style=node.right.incoming_edge_style())
    return node.pos

# scales the "outcomes" variables to be between 0 and 'new_range' for color selection
def scale_outcomes(outcomes, new_range=100):
    min_out = min(outcomes)
    max_out = max(outcomes)
    orig_range = round(max_out - min_out, 2)
    return [int(new_range * ((outcomes[i] - min_out) / orig_range))
        for i in range(len(outcomes))], min_out, max_out

# CALL THIS
# handles creating of the tree (the whole process)
def create_tree(numVerts, leafNames, homoVars, fname, outcomes=None, outcomeStr=None):
    # remove this if statement if you want an outcomeStr without outcomes
    if outcomes == None:
        outcomeStr = None
    # check input
    if len(numVerts) != len(leafNames):
        print 'Number of node values must match the number of leaves.'
        exit(-1)
    if len(homoVars) != len(leafNames):
        print 'Number of "homoVars" must match number of leaves'
        exit(-1)
    if outcomes != None and len(outcomes) != len(leafNames) + 1:
        print 'Number of outcomes must be one greater than the number of leaves'
        exit(-1)
    # scale outcome numbers to match colorbar
    min_out, max_out = None, None
    if outcomes != None:
        outcomes, min_out, max_out = scale_outcomes(outcomes)
    # map from leaf position to numeric value
    leaves_map = {leafNames[n]: numVerts[n] for n in range(len(leafNames))}
    # map from leaf position to homogenousness
    homo_map = {leafNames[n]: homoVars[n] for n in range(len(leafNames))}
    outcome_map = None
    if outcomes != None:
        # map from leaf position to outcome
        outcome_map = {leafNames[n]: outcomes[n + 1] for n in range(len(leafNames))}

    # reinitialize globals
    global cur_node_num
    global node_names
    cur_node_num = 0
    node_names = {}

    # build tree structure
    if outcomes != None:
        root = Node('g', outcome=outcomes[0])
    else:
        root = Node('g')
    for leaf in leafNames:
        while not in_tree(root, leaf):
            close = closest(root, leaf)
            nxt = direction(close, leaf)
            if nxt == 0:
                close.left = Node(close.pos + str(nxt), parent=close)
            else:
                close.right = Node(close.pos + str(nxt), parent=close)

    # set values
    set_values(root, leaves_map, outcome_map, homo_map)
    print_tree(root)

    # generate color spectrum values
    colors = [str(c) for c in list(Color("blue").range_to(Color("red"), 102))[1:-1]]
    # make colors valid length (unsimplified)
    for ndx in range(len(colors)):
        while len(colors[ndx]) < 7:
            colors[ndx] += '0'

    # build visualization
    graph = gv.Digraph(format=fname.split('.')[-1])
    draw_tree(graph, root, colors, leafNames)

    # render image
    graph = apply_styles(graph, STYLES)
    f = graph.render(filename='.'.join(fname.split('.')[:-1]))
    if outcomeStr != None or outcomes != None:
        # add gradient and title
        mti.modify(fname, outcomeStr, min_out, max_out)
    print '\nFinished creating', f

# if run from command line
if __name__ == '__main__':
    # read user input
    try:
        # leaf numerical values
        numVerts = [float(i.strip()) for i in sys.argv[1].split(',')]
        # leaf positions
        leafNames = [s.strip() for s in sys.argv[2].split(',')]
        # whether or not leaf defined by homogenous variable(s)
        homoVars = [int(i.strip()) == 1 for i in sys.argv[3].split(',')]
        # image filename
        fname = sys.argv[4]
        if len(sys.argv) > 5:
            # numbers that determine node color
            outcomes = [float(i.strip()) for i in sys.argv[5].split(',')]
            # text below color bar
            outcomeStr = sys.argv[6]
    except:
        print 'Usage: python create_tree_image.py "1,2,3,..." "g##,g###,g##,..." \
            "0,1,0,..." "0.9,.1,0.2,..." "this is the title" filename.png'
        exit(-1)
    # create tree
    create_tree(numVerts, leafNames, homoVars, fname, outcomes, outcomeStr)
