import matplotlib.pyplot as plt
import numpy

one_cns = []
two_cns = []
one_cns_upper = []
one_cns_lower = []
two_cns_upper = []
two_cns_lower = []

data = numpy.loadtxt(snakemake.input[0], usecols = (1, 2, 3))
x_p = snakemake.config["conta"]  

for i in range(0, len(data) - 1 , 2):
    one_cns.append((data[i][0]))
    two_cns.append((data[i+1][0]))
    one_cns_lower.append((data[i][1]))
    two_cns_lower.append((data[i+1][1]))
    one_cns_upper.append((data[i][2]))
    two_cns_upper.append((data[i+1][2]))

fig, ax = plt.subplots(1, 1)
plt.plot(x_p, one_cns, 'g', label = 'one-cns') 
for i in range(len(x_p)):
    line1 = plt.vlines(x_p[i], one_cns_lower[i], one_cns_upper[i])
    line1.set_color("green")
 
plt.plot(x_p, two_cns, 'b', label = 'two-cns') 
for i in range(len(x_p)):
    line1 = plt.vlines(x_p[i], two_cns_lower[i], two_cns_upper[i])
    line1.set_color("blue")
    
plt.plot(x_p, x_p, 'y--')

plt.legend()
plt.xlabel('c (expected)') 
plt.ylabel('c (estimated)')
plt.title('DoC = 1.00 x') 
plt.savefig(snakemake.output[0])
