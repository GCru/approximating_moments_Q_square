## Hall-Buckley-Eagleson method
import numpy
from momentchi2 import hbe


coeff_list=[0.4,0.3,0.3]
#coeff_list=[0.4,0.3,0.3]
t=2
coef_list=numpy.array(coeff_list)
x=hbe(coeff=coeff_list, x=t)
print(x)


coeff_list=[0.4,0.0001,0.0001]
t=2
coef_list=numpy.array(coeff_list)
x=hbe(coeff=coeff_list, x=t)
print(x)

exit()
# should give value close to 0.95, actually 0.94908
coeff_list=[0.8,0.1,0.1]
#coeff_list=[0.4,0.3,0.3]
t=4
coef_list=numpy.array(coeff_list)
x=hbe(coeff=coeff_list, x=t)

approx = t**(1/3) -(1-(2/9)*numpy.dot(coef_list, coef_list))
approx=approx/( ((2/9)*numpy.dot(coef_list, coef_list))**0.5)
print('*****', approx)
print(numpy.dot(coef_list, coef_list),numpy.dot(coeff_list,coeff_list)**0.5)

print(x)
exit()
coeff_list=[1/10]*10
x=hbe(coeff=coeff_list, x=2)
print(x)


#coeff_list=[0.005]*99
#coeff_list.append(1-0.005*99)
coeff_list=[0.5,0.1,0.1,0.1,0.1,0.1]
print(coeff_list, sum(coeff_list))
x=hbe(coeff=coeff_list, x=2)
print(x)

coeff_list=[0.5,0.11,0.09,0.1,0.1,0.1]
print(coeff_list, sum(coeff_list))
x=hbe(coeff=coeff_list, x=2)
print(x)
