using LinearAlgebra;

@enum Dof begin
	U
	P
	V
end

struct Node
	x::Float64;
	y::Float64;
	dofs::Set{Dof};
	vals::Dict{Dof, Float64};
end

function Node(x::Float64, y::Float64)
	dofs = Set{Dof}([U::Dof, V::Dof, P::Dof]);
	vals = Dict{Dof, Float64}(U::Dof => 0.0, V::Dof => 0.0, P::Dof => 0.0);
	return Node(x, y, dofs, vals);
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
	property::Property;
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
	nnode = getNodeNum(self);
	Ke = zeros(Float64, ndof, ndof);
	nu = self.property.nu;
	rho = self.property.rho;

	for p in self.points
		Der = getDerMatrix(self, p);
		N = getShapeMatrix(self, p);
		Nx = zeros(1, nnode);
		Ny = zeros(1, nnode);
		u0 = 0.0;
		v0 = 0.0;
		for i = 1 : nnode
			node = self.nodes[i];
			u0 += node.vals[U::Dof] * N[1, i];
			v0 += node.vals[V::Dof] * N[1, i];
			i += 1;
		end
		for i = 1 : nnode
			Nx[1, i] = Der[1, i];
			Ny[1, i] = Der[2, i];
		end
		
		J = getJacobi(self, p);
		k1 = 1;
		k2 = nnode + 1;
		k3 = nnode * 2 + 1;
		k4 = nnode * 3 + 1;
		w = abs(det(J)) * last(p);

		Ke[k1 : k2 - 1, k1 : k2 - 1] += (N' * u0 * Nx + N' * v0 * Ny + nu * Nx' * Nx + nu * Ny' * Ny) * w;
		Ke[k1 : k2 - 1, k2 : k3 - 1] += (N' * Nx / rho) * w;
		Ke[k2 : k3 - 1, k1 : k2 - 1] += (N' * Nx) * w;
		Ke[k2 : k3 - 1, k3 : k4 - 1] += (N' * Ny) * w;
		Ke[k3 : k4 - 1, k2 : k3 - 1] += (N' * Ny / rho) * w;
		Ke[k3 : k4 - 1, k3 : k4 - 1] += (N' * u0 * Nx + N' * v0 * Ny + nu * Nx' * Nx + nu * Ny' * Ny) * w;
	end
	return Ke;
end

function Tria3Element(nodes::Array{Node}, property::Property)
	return Tria3Element(nodes, [[1 / 3, 1 / 3, 1 / 3, 1 / 2]], property);
end

function getShapeMatrix(self::Tria3Element, p::Array{Float64})
	(l1, l2, l3, w) = p;
	nnode = getNodeNum(self);
	N = zeros(Float64, 1, nnode);
	for i = 1 : nnode
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
