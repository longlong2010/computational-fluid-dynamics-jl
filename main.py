from load import *;
from geometry import *;
from model import *;
import sys
sys.path.append('pyNastran');
from pyNastran.bdf.bdf import BDF, read_bdf;
if __name__ == '__main__':
	numpy.set_printoptions(threshold = numpy.nan);

	model = Model();
	bdf = BDF(debug=False);
	bdf.read_bdf('A1.bdf', punch=True);

	for x in bdf.spcs:
		spc = bdf.get_spcs(x);


	nodes = dict();
	for x in bdf.nodes:
		n = bdf.nodes[x];
		c = n.xyz;
		node = Node(c[0], c[1]);
		if x in spc[0]:
			node.addConstraint(Constraint.PHI);
		if c[0] < 1e-5 or abs(c[0] - 5) < 1e-5:
			node.addConstraint(Constraint.PHI, c[1] - 0.5);
		if c[1] < 1e-5:
			node.addConstraint(Constraint.PHI, -0.5);
		if abs(c[1] - 1) < 1e-5:
			node.addConstraint(Constraint.PHI, 0.5);
		nodes[x] = node;

	for x in bdf.elements:
		e = bdf.elements[x];
		[i, j, k] = e.nodes;
		n1 = nodes[i];
		n2 = nodes[j];
		n3 = nodes[k];
		
		model.addElement(Tria3Element(n1, n2, n3));
	
	model.solve();
