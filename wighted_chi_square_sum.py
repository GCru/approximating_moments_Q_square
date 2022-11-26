## Hall-Buckley-Eagleson method
from momentchi2 import hbe

# should give value close to 0.95, actually 0.94908
x=hbe(coeff=[1.5, 1.5, 0.5, 0.5], x=10.203)

print(x)

coeff_list=[1/10]*10
x=hbe(coeff=coeff_list, x=2)
print(x)


coeff_list=[0.005]*99
coeff_list.append(1-0.005*99)
print(coeff_list)
x=hbe(coeff=coeff_list, x=2.5)
print(x)


# should give value close to 0.95, actually 0.94908
x=hbe(coeff=[0.25,0.25], x=1.0)

print(x)