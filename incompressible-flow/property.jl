abstract type Property end;

struct Material <: Property
	nu::Float64;
	rho::Float64;
end
