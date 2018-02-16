import abc;
import numpy;
from enum import Enum;
Dof = Enum('Dof', 'X Y PHI');
from load import *;

class Node:
	def __init__(self, x, y):
		self.x = x;
		self.y = y;

		self.dofs = set();
		self.dofs.add(Dof.PHI);
		
		self.loads = dict();
		self.constraints = dict();
		self.values = dict();

	def addLoad(self, load, v):
		if self.loads.__contains__(load):
			self.loads[load] += v;
		else:
			self.loads[load] = v;
		
	def addConstraint(self, constraint, v = 0):
		self.constraints[constraint] = v;
	
	def getLoads(self):
		return self.loads;

	def setValue(self, dof, v):
		if dof in self.dofs:
			self.values[dof] = v;

	def getConstraints(self):
		return self.constraints;

	def getDofs(self):
		dofs = [];
		for dof in Dof:
			if dof in self.dofs:
				dofs.append(dof);
		return dofs;

	def getDofNum(self):
		return len(self.dofs);

class Element(metaclass = abc.ABCMeta):
	def __init__(self, nodes, points):
		self.nodes = nodes;
		self.points = points;
	
	def getNodeNum(self):
		return len(self.nodes);

	def getDofNum(self):
		ndof = 0;
		for node in self.nodes:
			ndof += node.getDofNum();
		return ndof;
	
	def getNodes(self):
		return self.nodes;

	def getCoord(self):
		nnode = self.getNodeNum();
		coord = numpy.zeros((nnode, 2));
		i = 0;
		for n in self.nodes:
			coord[i] = numpy.array([n.x, n.y]);
			i += 1;
		return coord;

	@abc.abstractmethod	
	def getShapeDerMatrix(self, p):
		pass;
	
	@abc.abstractmethod
	def getShapeMatrix(self, p):
		pass;
	
	def getJacobi(self, p):
		der = self.getShapeDerMatrix(p);
		coord = self.getCoord();
		return der.dot(coord);

	def getVelocityMatrix(self, p):
		der = self.getShapeDerMatrix(p);
		J = self.getJacobi(p);
		Der = numpy.linalg.inv(J).dot(der);
		return Der;

	def getStiffMatrix(self):
		ndof = self.getDofNum();
		Ke = numpy.zeros((ndof, ndof));
		for p in self.points:
			B = self.getVelocityMatrix(p);
			J = self.getJacobi(p);
			Ke += B.T.dot(B) * abs(numpy.linalg.det(J)) * p[3];
		return Ke;
	
class Tria3Element(Element):
	def __init__(self, n1, n2, n3):
		super(Tria3Element, self).__init__([n1, n2, n3], [[1 / 3, 1 / 3, 1 / 3, 1 / 2]]);

	def getShapeMatrix(self, p):
		[l1, l2, l3, w] = p;
		ndof = self.getDofNum();
		N = numpy.array([[l1, l2, l3]]);
		return N;

	def getShapeDerMatrix(self, p):
		return numpy.array([
			[1, 0, -1],
			[0, 1, -1],
		]);
