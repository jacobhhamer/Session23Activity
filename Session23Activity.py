import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import astropy.units as unit
import astropy.constants as cons
import matplotlib as mpl
mpl.rcParams['font.size']=16

def f(r):
    x,y,xdot,ydot=r[2],r[3],x2(r[2],r[3]),y2(r[2], r[3])
    return np.array([x,y,xdot,ydot])

def trajectory(v0, theta, planet, mass=1.):
    if (planet=='Earth' or planet=='earth'):
        g=9.8
        rho=1.22
    elif (planet=='Mars' or planet=='mars'):
        g=3.71
        rho=0.20
    else: 
        print('Planets currently supported are Earth and Mars. Please enter one of those two for planet.')
        return
    
    m=1.
    R=.08
    theta=np.deg2rad(theta)
    x0,x1,y0,y1=0.,v0*np.cos(theta), 0, v0*np.sin(theta)
    xs=[]
    ys=[]
    C=0.47
    H=.01
    t0=0
    delta=1e-4

    y2=lambda x1,y1: -g-(np.pi*R**2*rho*C*y1*np.sqrt(x1**2+y1**2))/(2*m)
    x2=lambda x1,y1: -(np.pi*R**2*rho*C*x1*np.sqrt(x1**2+y1**2))/(2*m)



    vec=np.array([x0,y0,x1,y1])
    y=y0
    while y>-1e-6:
        xs.append(vec[0])
        ys.append(vec[1])
        n=1

        r1=vec+0.5*H*f(vec)
        r2=vec+H*f(r1)

        R1=np.empty([1,4], float)
        R1[0]=0.5*(r1+r2+0.5*H*f(r2))

        error=2*H*delta
        while error>H*delta:
            n+=1
            h=H/n

            r1=vec+0.5*h*f(vec)
            r2=vec+h*f(r1)   
            for i in range(n-1):
                r1+=h*f(r2)
                r2+=h*f(r1)

            R2=R1
            R1=np.empty([n,4],float)
            R1[0]=0.5*(r1+r2+0.5*h*f(r2))

            for m in range(1,n):
                epsilon=(R1[m-1]-R2[m-1])/((n/(n-1))**(2*m)-1.)
                R1[m]=R1[m-1]+epsilon
            error=np.abs(epsilon[0])
        vec=R1[n-1]
        y=ys[-1]
    return xs, ys

 def max_distance(v0, theta, planet, mass=1., plot=False):
 	xs, ys = trajectory(v0=v0, theta=theta, planet=planet)
 	if plot:
	 	f, ax = plt.subplots(1, figsize=(7,5))
	 	ax.plot(xs, ys)
	 	ax.set_xlabel('x[m]')
	 	ax.set_ylabel('y[m]')
	 	ax.set_title('$v_0={}$, $\\theta={}$, {}, $x_{max}=${:.2f}'.format(v0,theta, planet, xs[-1]))
 		plt.show()
 	return xs[-1]

 def multimass_distance(v0, theta, planet, minmass, maxmass):
 	masses=np.linspace(minmass, maxmass, 100)
 	xlist=[]
 	for m in masses:
 		x=max_distance(v0=v0, theta=theta, planet=planet, mass=m)
 		xlist.append(x)
 	return masses, xlist

 def multimass_distance_plot(v0, theta, planet, minmass, maxmass):
 	masses, xlist=multimass_distance(v0=v0, theta=theta, planet=planet, minmass=minmass, maxmass=maxmass)
 	f, ax = plt.subplots(1, figsize=(7,5))
 	ax.plot(masses, xlist)
 	ax.set_ylabel('Mass [kg]')
 	ax.set_xlabel('Maximum Distance [m]')
 	plt.show()
 	return f, ax