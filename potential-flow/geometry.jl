using LinearAlgebra;

@enum Dof begin
	PHI
end

struct Node
	x::Float64;
	y::Float64;
	dofs::Set{Dof};
end

function Node(x::Float64, y::Float64)
	dofs = Set{Dof}([PHI::Dof]);
	return Node(x, y, dofs);
end

function getDofNum(self::Node)
	return length(self.dofs);
end

abstract type Element end;
abstract type Element2D <: Element end;
abstract type TriaElement <: Element2D end;

mutable struct Tria3Element <: TriaElement
	nodes::Array{Node};
	points::Array{Array{Float64}};
end

function getNodeNum(self::Element)
	length(self.nodes);
end

function getDofNum(self::Element)
	ndof::Int32 = 0;
	for node in self.nodes
		ndof += getDofNum(node);
	end
	return ndof;
end

function getCoord(self::Element2D)
	nnode = getNodeNum(self);
	coord = zeros(Float64, nnode, 2);
	i::Int32 = 1;
	for node in self.nodes
		coord[i, :] = Array{Float64}([node.x, node.y]);
		i += 1;
	end
	return coord;
end

function getDerMatrix(self::Element2D, p::Array{Float64})
	J = getJacobi(self, p);
	der = getShapeDerMatrix(self, p);
	return J^-1 * der;
end;

function getJacobi(self::Element2D, p::Array{Float64})
	der = getShapeDerMatrix(self, p);
	coord = getCoord(self);
	return der * coord;
end

function getStiffMatrix(self::Element2D)
	ndof = getDofNum(self);
	Ke = zeros(Float64, ndof, ndof);
	for p in self.points
		B = getDerMatrix(self, p);
		J = getJacobi(self, p);
		Ke += B' * B * abs(det(J)) * last(p);
	end
	return Ke;
end

function Tria3Element(nodes::Array{Node})
	return Tria3Element(nodes, [[1 / 3, 1 / 3, 1 / 3, 1 / 2]]);
end

function getShapeMatrix(self::Tria3Element, p::Array{Float64})
	(l1, l2, l3, w) = p;
	ndof = getDofNum(self);
	N = zeros(Float64, 1, ndof);
	for i = 1 : ndof
		N[1, i] = p[i];
	end
	return N;
end

function getShapeDerMatrix(self::Tria3Element, p::Array{Float64})
	nnode = getNodeNum(self);
	der = zeros(Float64, 2, nnode);
	der[1, 1] = 1;
	der[2, 2] = 1;
	der[1, 3] = -1;
	der[2, 3] = -1;
	return der;
end
