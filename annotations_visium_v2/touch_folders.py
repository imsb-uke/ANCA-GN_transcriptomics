import os

vs = range(1,7)
rs = ["A", "B", "C", "D"]

for v in vs:
	for r in rs:
		folder = f"V{v}_{r}"
		if not os.path.exists(folder):
			os.mkdir(folder)

