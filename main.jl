include("geometry.jl");
include("load.jl");
include("model.jl");
begin
	local n1 = Node(0.0, 1.0);
	local n2 = Node(0.0, 0.0);
	local n3 = Node(1.0, 0.0);
	local n4 = Node(1.0, 1.0);
	local n5 = Node(2.0, 0.0);
	local n6 = Node(2.0, 1.0);

	local e1 = Tria3Element([n1, n2, n3]) 
	local e2 = Tria3Element([n1, n3, n4]);
	local e3 = Tria3Element([n3, n4, n5]);
	local e4 = Tria3Element([n4, n5, n6]);

	local spc = SPC();
	addConstraint(spc, PHI::Dof);
	addNode(spc, n5);
	addNode(spc, n6);

	local load = Load(0.5);
	addNode(load, n1);
	addNode(load, n2);
	
	local m = Model();
	addElement(m, e1);
	addElement(m, e2);
	addElement(m, e3);
	addElement(m, e4);
	addConstraint(m, spc);
	addLoad(m, load);
	solve(m);

	local nodes::Dict{Int32, Node} = Dict{Int32, Node}();
	local elements::Array{Array{Int32}} = Array{Array{Int32}}([]);
	local loads::Dict{Int32, Float64} = Dict{Int32, Float64}();
	
	m = Model();

	for line in eachline("A7.bdf")
		if line[1] != '$'
			local len = length(line);
			local n = len ÷ 8 + 1;
			local card = strip(line[1 : 8]);
			if card == "CTRIA3"
				local nids::Array{Int32} = Array{Int32}([]);
				local id = parse(Int32, String(strip(line[9 : 16])));
				for nid in split(strip(line[25 : end]), r"\s+")
					push!(nids, parse(Int32, String(nid)));
				end
				push!(elements, nids);
			elseif card == "GRID"
				local x::Float64 = parse(Float64, replace(String(strip(line[25 : 32])), r"\-(\d+)$" => s"e-\1"));
				local y::Float64 = parse(Float64, replace(String(strip(line[33 : 40])), r"\-(\d+)$" => s"e-\1"));
				local id = parse(Int32, String(strip(line[9 : 16])));
				node = Node(x, y);
				nodes[id] = node;
			end
		end
	end

	local l1 = Load(0.1);
	local l2 = Load(0.05);
	for i = 21 : 29
		addNode(l1, nodes[i]);
	end
	addNode(l2, nodes[2]);
	addNode(l2, nodes[20]);
	
	spc = SPC();
	addConstraint(spc, PHI::Dof);
	addNode(spc, nodes[1]);
	addNode(spc, nodes[30]);
	for i = 48 : 56
		addNode(spc, nodes[i]);
	end

	for nids in elements
		local (i, j, k) = nids;
		local n1::Node = nodes[i];
		local n2::Node = nodes[j];
		local n3::Node = nodes[k];

		local e = Tria3Element([n1, n2, n3]);
		addElement(m, e);		
	end
	addLoad(m, l1);
	addLoad(m, l2);
	addConstraint(m, spc);
	solve(m);
	save(m, "data");
end
