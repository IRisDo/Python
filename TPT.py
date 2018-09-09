import numpy as np
import matplotlib.pyplot as plt

#variable initial
phi1 = 1
W0 = 0.000317
S = 0.002376
ld = 0.543
m1 = 0.00039
m2 = 0.00137
p0 = 15000
om = 0.00008
f = 950000
alpha = 1
delta = 1.6
Ik = 62.5
kh = 1.62
teta = 0.12
h = 0.000001
muy = 0.06
lamda = -0.44
g = 98.1
s1 = 1
s2 = 1
phi=phi1+om/((m1+m2)*3)
p = p0
pmax = 0

#create value for solve
n = 4
pt = np.zeros((1,4))    #create array zeros(1,4)
x = np.zeros((1,5))     #create array zeros(1,5)
y = np.zeros((1,5))      #create array zeros(1,5)

#Set value intial
def set_value_intial():
	psi0=(W0/om-1/delta)/(f/p0+alpha-1/delta)
	def found_z0():
		hs = [kh*muy, kh*lamda, kh, -psi0]
		gt = np.roots(hs)
		for i in range(0,len(gt)):
			if (gt[i].imag == 0) and (gt[i].real <= 1) and (gt[i].real > 0):
				return gt[i].real

	x[0][0] = found_z0()     #z
	x[0][1] = psi0           #psi
	x[0][2] = 0              #v
	x[0][3] = 0              #l
	x[0][4] = 0              #t
	

def hptvp():
	pt[0][0] = s1*p/Ik
	if muy == 0:
		pt[0][1] = sqrt((kh**2)-4*(kh-1)*x[0][1])*pt[0][0]
	else:
		pt[0][1] = kh*(1+x[0][0]*(2*lamda+3*muy*x[0][0]))*pt[0][0]
	pt[0][2] = s2*S*p*g/(phi*(m1+m2))
	pt[0][3] = x[0][2]

def rkut():
	k = np.zeros((n,4))        #create array zeros(n,4)
	def solve_k():
		hptvp()
		for i in range(0,n):
			k[i][s] = pt[0][i]*h

	for i in range(0,n+1):
		y[0][i] = x[0][i]
	for s in range(0,4):
		if s == 0:
			solve_k()
			for i in range(0,n):
				x[0][i] = y[0][i] + k[i][s]/2
		elif s == 1:
			solve_k()
			for i in range(0,n):
				x[0][i] = y[0][i] + k[i][s]/2
		elif s == 2:
			solve_k()
			for i in range(0,n):
				x[0][i] = y[0][i] + k[i][s]
		else:
			solve_k()
			for i in range(0,n):
				x[0][i] = y[0][i] + (k[i][0] + k[i][3])/6 + (k[i][1]+k[i][2])/3
	x[0][n] += h

def raw():
	# subplot(121)
	fig, ax1 = plt.subplots()

	color = 'tab:red'
	ax1.set_xlabel('time(ms)')
	ax1.set_ylabel('p(KG/cm^2)', color=color)
	plt.title('Do thi ap suat va van toc theo thoi gian')
	ax1.plot(KQ[:,4], KQ[:,5], color=color)
	ax1.tick_params(axis='y', labelcolor=color)

	ax2 = ax1.twinx()

	color = 'tab:blue'
	ax2.set_ylabel('v(m/s)', color=color)
	ax2.plot(KQ[:,4], KQ[:,2],color=color)
	ax2.tick_params(axis='y', labelcolor=color)
	fig.tight_layout()

	# subplot(122)
	fig, ax3 = plt.subplots()
	color = 'tab:red'
	ax3.set_xlabel('length(dm)')
	ax3.set_ylabel('p(KG/cm^2)', color=color)
	plt.title('Do thi ap suat va van toc theo quang duong')
	ax3.plot(KQ[:,3], KQ[:,5], color=color)
	ax3.tick_params(axis='y', labelcolor=color)
	
	ax4 = ax3.twinx()

	color = 'tab:blue'
	ax4.set_ylabel('v(m/s)', color=color)
	ax4.plot(KQ[:,3], KQ[:,2],color=color)
	ax4.tick_params(axis='y', labelcolor=color)
	fig.tight_layout()

	plt.show()

if __name__ == "__main__":
	set_value_intial()
	print(x[0][1])
	KQ = np.array([[x[0][0], x[0][1], x[0][2]/10, x[0][3], x[0][4]*1000, p/100]])       #list contain solution
	while x[0][3] < ld:
		if x[0][0] < 1:
			s1 = 1;
			tk = x[0][4]*1000
			pk = p/100
			lk = x[0][3]
			vk = x[0][2]/10
		else:
			s1 = 0
			x[0][0] = 1
			x[0][1] = kh*(1+lamda+muy)
		if pmax < p:
			pmax = p
			tm = x[0][4]*1000
			lm = x[0][3]
			vm = x[0][2]/10
		rkut()
		p = (f*om*x[0][1]-teta*phi*(m1+m2)*(x[0][2]**2)/(2*g))/(W0-om*(1-x[0][1])/delta-alpha*om*x[0][1]+S*x[0][3])
		KQ = np.append(KQ, [[x[0][0], x[0][1], x[0][2]/10, x[0][3], x[0][4]*1000, p/100]], axis=0)
	print(KQ)
	raw()
