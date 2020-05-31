from matplotlib.tri import Triangulation, CubicTriInterpolator;
import matplotlib.pyplot as plt;
import numpy;


if __name__ == '__main__':
	x = [];
	y = [];
	v = [];
	m = [];
	for line in open('data.1'):
		row = line.strip().split(" ");
		x.append(float(row[0]));
		y.append(float(row[1]));
	for line in open('data.2'):
		row = line.strip().split(" ");
		v.append(float(row[0]));
	for line in open('data.3'):
		row = line.strip().split(" ");
		m.append([int(row[0]) - 1, int(row[1]) - 1, int(row[2]) - 1]);
	
	triangles = Triangulation(x, y, m);

	tci = CubicTriInterpolator(triangles, -numpy.array(v));
	(u, v) = tci.gradient(triangles.x, triangles.y);
	v_norm = numpy.sqrt(u ** 2, v ** 2);

	fig, ax = plt.subplots();
	ax.set_aspect('equal');
	ax.use_sticky_edges = False;
	ax.margins(0.07);
	ax.triplot(triangles, color='0.8', lw=0.5);

	#tcf = ax.tricontourf(triangles, -numpy.array(v));
	#fig.colorbar(tcf);

	ax.quiver(triangles.x, triangles.y, u / v_norm, v / v_norm, color='blue');
	plt.savefig('1.svg');
