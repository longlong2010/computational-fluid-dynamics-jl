from geometry import *;
import scipy.sparse;

class Model:
	def __init__(self):
		self.nodes = [];
		self.elements = [];

	def addElement(self, e):
		nodes = e.getNodes();
		for n in nodes:
			if n not in self.nodes:
				self.nodes.append(n);
		if e not in self.elements:
			self.elements.append(e);

	def getDofNum(self):
		ndof = 0;
		for node in self.nodes:
			ndof += node.getDofNum();
		return ndof;

	def getElements(self):
		return self.elements;

	def init(self):
		ndof = self.getDofNum();
		self.K = numpy.zeros((ndof, ndof));
		self.R = numpy.zeros((ndof, 1));

	def integrate(self):
		ndof = self.getDofNum();
		for e in self.elements:
			nodes = e.getNodes();
			size = e.getDofNum();
		
			m = dict();
			l = 0;
			for n in nodes:
				k = self.nodes.index(n) * n.getDofNum();
				for i in range(0, n.getDofNum()):
					m[l] = k + i;
					l += 1;
			Ke = e.getStiffMatrix();
			for i in range(0, size):
				for j in range(0, size):
					mi = m[i];
					mj = m[j];
					self.K[mi][mj] += Ke[i][j];

	def integrateLoad(self):
		k = 0;
		for n in self.nodes:
			dofs = n.getDofs();
			dofn = 0;
			for d in dofs:
				loads = n.getLoads();
				for l, v in loads.items():
					if l.getDof() == d:
						self.R[k + dofn][0] = v;
				dofn += 1;
			k += n.getDofNum();


	def addConstraint(self):
		ndof = self.getDofNum();
		k = 0;
		for n in self.nodes:
			constraints = n.getConstraints();
			dofs = n.getDofs();

			dofn = 0;
			for d in dofs:
				for c, v in constraints.items():
					if c.getDof() == d:
						self.K[k + dofn][k + dofn] += 1e20;
						self.R[k + dofn][0] = v * self.K[k + dofn][k + dofn];
				dofn += 1;
			k += n.getDofNum();
	
	def solveEquations(self):
		u = numpy.linalg.solve(self.K, self.R);
		k = 0;
		for n in self.nodes:
			dofs = n.getDofs();
			dofn = 0;
			for d in dofs:
				n.setValue(d, u[k + dofn][0]);
				dofn += 1;
			k += n.getDofNum();

	def outputResult(self):
		for e in self.elements:
			stress = e.getStress();
			print(stress[0][0]);

	def solve(self):
		self.init();
		self.integrate();
		self.integrateLoad();
		self.addConstraint();
		self.solveEquations();
		self.outputResult();	
