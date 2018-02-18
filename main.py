from load import *;
from geometry import *;
from model import *;
from pyNastran.bdf.bdf import BDF, read_bdf;
if __name__ == '__main__':
	numpy.set_printoptions(threshold = numpy.nan);

	model = Model();
	bdf = BDF(debug=False);
	bdf.read_bdf('A1.bdf', punch=True);

	for x in bdf.spcs:
		spc = bdf.get_spcs(x);

	for x in bdf.elements:
		e = bdf.elements[x];
		[i, j, k] = e.nodes;
		c1 = bdf.nodes[i].xyz;
		c2 = bdf.nodes[j].xyz;
		c3 = bdf.nodes[k].xyz;

		n1 = Node(c1[0], c1[1]);
		n2 = Node(c2[0], c2[1]);
		n3 = Node(c3[0], c3[1]);
		if i in spc[0]:
			n1.addConstraint(Constraint.PHI);
		if j in spc[0]:
			n2.addConstraint(Constraint.PHI);
		if k in spc[0]:
			n3.addConstraint(Constraint.PHI);
		model.addElement(Tria3Element(n1, n2, n3));
	
	model.solve();
