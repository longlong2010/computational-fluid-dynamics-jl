include("property.jl");
include("geometry.jl");
include("load.jl");
include("model.jl");
begin
	local nodes::Dict{Int32, Node} = Dict{Int32, Node}();
	local elements::Array{Array{Int32}} = Array{Array{Int32}}([]);
	local loads::Dict{Int32, Float64} = Dict{Int32, Float64}();
	
	local mat = Material(0.1, 1.0);

	local m = Model();

	for line in eachline("A9.bdf")
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


	local spc1 = SPC();
	addConstraint(spc1, U::Dof, 5.0);
	addConstraint(spc1, V::Dof, 0.0);
	addConstraint(spc1, P::Dof, 0.0);
	for i = 92 : 100
		addNode(spc1, nodes[i]);
		nodes[i].vals[U::Dof] = 5.0;
	end
	for i = 1 : 91
		addNode(spc1, nodes[i]);
		nodes[i].vals[U::Dof] = 5.0;
	end

	local spc2 = SPC();
	addConstraint(spc2, U::Dof, 0.0);
	addConstraint(spc2, V::Dof, 0.0);
	addConstraint(spc2, P::Dof, 0.0);
	for i = 101 : 150
		addNode(spc2, nodes[i]);
	end

	for nids in elements
		local (i, j, k) = nids;
		local n1::Node = nodes[i];
		local n2::Node = nodes[j];
		local n3::Node = nodes[k];

		local e = Tria3Element([n1, n2, n3], mat);
		addElement(m, e);		
	end
	addConstraint(m, spc1);
	addConstraint(m, spc2);
	local n = parse(Int32, ARGS[1]);
	for i = 1 : n
		solve(m, 0.05);
		save(m, "data");
		run(`python main.py $i`);
		if i < 10
			run(`mv $i.png 00$i.png`);
		elseif i < 100
			run(`mv $i.png 0$i.png`);
		end
	end
	run(`convert -delay 10 \*.png -loop 0 animated.gif`);
end
