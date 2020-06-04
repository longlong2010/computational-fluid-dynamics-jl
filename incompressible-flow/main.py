from matplotlib.tri import Triangulation, CubicTriInterpolator;
import matplotlib.pyplot as plt;
import numpy;


if __name__ == '__main__':
	x = [];
	y = [];
	u = [];
	v = [];
	m = [];
	for line in open('data.1'):
		row = line.strip().split(" ");
		x.append(float(row[0]));
		y.append(float(row[1]));
	for line in open('data.2'):
		row = line.strip().split(" ");
		u.append(float(row[0]));
		v.append(float(row[2]));
	for line in open('data.3'):
		row = line.strip().split(" ");
		m.append([int(row[0]) - 1, int(row[1]) - 1, int(row[2]) - 1]);

	fig, ax = plt.subplots();
	ax.set_aspect('equal');
	q = ax.quiver(x, y, u, v);
	plt.savefig('1.svg');
