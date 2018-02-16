from load import *;
from geometry import *;
from model import *;

import csv;

if __name__ == '__main__':
	numpy.set_printoptions(threshold = numpy.nan);

	model = Model();
	

	nodes = dict();
	l = 0;
	with open('demo.msh', 'r') as f:
		reader = csv.reader(f, delimiter = ' ');
		mode = '';
		for row in reader:
			if row[0] == '$Nodes':
				mode = 'Node';
				continue;
			elif row[0] == '$Elements':
				mode = 'Element';
				continue;

			if mode == 'Node':
				if len(row) == 4:
					n = Node(float(row[1]) / 1000, float(row[2]) / 1000, float(row[3]) / 1000);
					if abs(n.x + 5e-2) < 1e-10:
						n.addConstraint(Constraint.X);
						n.addConstraint(Constraint.Y);
						n.addConstraint(Constraint.Z);
					elif abs(n.x - 5e-2) < 1e-10:
						n.addLoad(Load.X, 10);
						l += 10;
					nodes[row[0]] = n;
			elif mode == 'Element':
				if len(row) == 9:
					n1 = nodes[row[5]];
					n2 = nodes[row[6]];
					n3 = nodes[row[7]];
					n4 = nodes[row[8]];
					model.addElement(Tria3Element(n1, n2, n3));
	model.solve();
