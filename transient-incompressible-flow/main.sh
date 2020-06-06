for i in `seq 1 20`
do
	julia main.jl "$i"
	python main.py "$i"
done

for i in `seq 1 9`
do
	mv "$i.png" "0$i.png"
done
convert -delay 100 *.png -loop 0 animated.gif
