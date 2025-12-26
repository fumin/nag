numElements = 31
toN = x -> (lift(x, ZZ) + numElements) % numElements
toS = x -> toString(toN(x))

k = GF numElements
a0 = a - toN(a)
for i from 0 to numElements-1 do(
	ei = a0 + i;
	s = "{";
	for j from 0 to numElements-1 do(
		ej = a0 + j;

		comma = ", ";
		if j == 0 then comma = "";

		eij = ei * ej;
        	s = s | comma | toS(eij);
	);
	s = s | "},";
	print(s);
)
