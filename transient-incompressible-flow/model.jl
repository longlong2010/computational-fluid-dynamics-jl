mutable struct Model
	nodes::Array{Node};
	elements::Array{Element};
	constraints::Array{Constraint};
	result::Array{Float64};
end

function addElement(self::Model, e::Element)
	for node in e.nodes
		if !(node in self.nodes)
			push!(self.nodes, node);
		end
	end
	if !(e in self.elements)
		push!(self.elements, e);
	end
end

function addConstraint(self::Model, c::Constraint)
	push!(self.constraints, c);
end

function addLoad(self::Model, l::Load)
	push!(self.loads, l);
end

function getDofNum(self::Model)
	ndof::Int32 = 0;
	for node in self.nodes
		ndof += getDofNum(node);
	end
	return ndof;
end

function getNodeNum(self::Model)
	return length(self.nodes);
end

function solve(self::Model, dt::Float64)
	ndof = getDofNum(self);
	nnode = getNodeNum(self);
	K = zeros(Float64, ndof, ndof);
	M = zeros(Float64, ndof, ndof);
	R = zeros(Float64, ndof, 1);
	inode::Int32 = 1;
	mnode::Dict{Node, Int32} = Dict{Node, Int32}();
	for node in self.nodes
		mnode[node] = inode;
		R[inode] = node.vals[U::Dof];
		R[inode + nnode] = node.vals[P::Dof];
		R[inode + nnode * 2] = node.vals[V::Dof];
		inode += 1;
	end
	
	for element in self.elements
		minode::Dict{Int32, Int32} = Dict{Int32, Int32}();
		einode::Int32 = 1;
		ennode = length(element.nodes);
		for node in element.nodes
			inode = mnode[node];
			for dofn = 1 : getDofNum(node)
				minode[einode + (dofn - 1) * ennode] = inode + (dofn - 1) * nnode;
			end
			einode += 1;
		end

		Ke = getStiffMatrix(element);
		Me = getMassMatrix(element);
		(m, n) = size(Ke);
		for i = 1 : m
			mi = minode[i];
			for j = 1 : n
				mj = minode[j];
				K[mi, mj] += Ke[i, j];
				M[mi, mj] += Me[i, j];
			end
		end
	end

	R = (M - K * dt / 2) * R;
	K = (M + K * dt / 2);
	
	for constraint in self.constraints
		for node in constraint.nodes
			inode = mnode[node];
			dofn::Int32 = 1;
			for d in instances(Dof)
				if haskey(constraint.values, d)
					idof = inode + (dofn - 1) * nnode;
					K[idof, idof] += 1e20;
					R[idof] = constraint.values[d] * K[idof, idof];
				end
				dofn += 1;
			end
		end
	end

	self.result = K \ R;
	inode = 1;
	for node in self.nodes
		node.vals[U::Dof] = self.result[inode];
		node.vals[P::Dof] = self.result[inode + nnode];
		node.vals[V::Dof] = self.result[inode + nnode * 2];

		inode += 1;
	end
end

function save(self::Model, file::String)
	local io1 = open(file * ".1", "w");
	local io2 = open(file * ".2", "w");
	local io3 = open(file * ".3", "w");

	nnode::Int32 = getNodeNum(self);
	inode::Int32 = 1;
	for node in self.nodes
		local ndof = getDofNum(node);
		v = zeros(Float64, ndof, 1);
		for i = 1 : ndof
			idof = inode + (i - 1) * nnode;
			v[i] = self.result[idof];
		end
		inode += 1;
		write(io1, join([node.x, node.y], " ") * "\n");
		write(io2, join(v, " ") * "\n");
	end

	for element in self.elements
		local nids::Array{Int32} = Array{Int32}([]);
		for node in element.nodes
			local id::Int32 = findfirst(isequal(node), self.nodes);
			push!(nids, id);
		end
		write(io3, join(nids, " ") * "\n");
	end
	close(io1);
	close(io2);
	close(io3);
end

function Model()
	return Model(Array{Node}([]), Array{Element}([]), Array{Constraint}([]), Array{Float64}([]));
end
