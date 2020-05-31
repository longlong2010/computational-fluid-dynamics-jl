abstract type Constraint end;

struct SPC <: Constraint
	values::Dict{Dof, Float64};
	nodes::Set{Node};
end

struct Load
	values::Dict{Dof, Float64};
	nodes::Set{Node};
end

function addConstraint(self::SPC, dof::Dof)
	addConstraint(self, dof, 0.0);
end

function addConstraint(self::SPC, dof::Dof, val::Float64)
	self.values[dof] = val;
end

function addNode(self::SPC, node::Node)
	push!(self.nodes, node);
end

function SPC()
	return SPC(Dict{Dof, Float64}(), Set{Node}([]));
end

function addNode(self::Load, node::Node)
	push!(self.nodes, node);
end

function Load(lphi::Float64)
	return Load(Dict{Dof, Float64}(PHI::Dof => lphi), Set{Node}([]));
end
