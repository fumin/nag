-- M2 --script macaulay.m2

loadPackage "AssociativeAlgebras"
A = QQ<|a, b|>
I = ideal {
	a^2-1,
	b^3-1,
	(a*b*a*b*a*(b^2)*a*(b^2))^2-1
}

prevNumBases = -1
for i from 0 to 999999 do(
	Igb = NCGB(I, 999999);
	numBases = numgens(source(Igb));
	if numBases == prevNumBases then break;

	I = ideal(Igb);
	prevNumBases = numBases;
	print(i, numBases);
)

"basis.m2" << toString Igb << close
