numElements = 5
toN = x -> (lift(x, ZZ) + numElements) % numElements

k = GF numElements
a0 = a - toN(a)
R = k[x]

f = (x^2+1)*(x-2)
g = (x^3+7)*(x-2)
print("f", toString(f))
print("g", toString(g))
h = gcd(f, g)
print("gcd", toString(h))
