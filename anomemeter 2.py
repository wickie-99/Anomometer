import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def property(T):
        Tf = (T+300)/2
        density = 1.1614 + (Tf-300)/(350-300)*(0.995-1.1614)
        viscosity = (184.6 + (Tf-300)/(350-300)*(208.2-184.6))*10**-7
        conductivity = (26.3 + (Tf-300)/(350-300)*(30.0-26.3))*10**-3
        Pr = (0.707 + (Tf-300)/(350-300)*(0.70-0.707))
        return(density,viscosity,conductivity,Pr)

def Nu_cal(Re_num, Pr):

    if (Re_num<5*10**5):
        return 0.664*Re_num**0.5*Pr*0.33
    else:
        return (0.037*Re_num**0.8*-871)*Pr*0.33

I, k,Pr,s,l,Ts, Too = sp.symbols("I, k,Pr,s,l,T_s, T_oo")
Re, rho,t, mu, v, Nu = sp.symbols("Re, rho,t, mu, v, Nu")
dI, dv = sp.symbols("dI, dv")
Q1, Q2 = sp.symbols("Q_1, Q_2")
Err = sp.Symbol("Err")
w = sp.Symbol('w')

eq_Re = sp.Eq(Re,rho*v*w/mu)
eq_Q1 = sp.Eq(Q1, I**2*s*l/w/t)
eq_Q2 = sp.Eq(Q2, k/w*Nu *2*w*l*(Ts-Too))

l_num = 0.01
s_num = 10.6 * 10**-8
Too_num = 300

eq_Q1
v_max = 50
I_max = 0.2
w_max = 2
temp_acc = 0.1
w_acc = 0.1

values_1 = {v:v_max, I:I_max, l:l_num,s:s_num, Too:Too_num}

Error = 10**10
W = 0
temperature = 0

T_step = int((327-312)/temp_acc)
w_step = int(2/w_acc)


for temp in np.linspace(313,327,T_step):
    prop = property(temp)
    values_2 = {k:prop[2], mu:prop[1], rho: prop[0], Pr:prop[3]}
    for width in np.linspace(0.1,2,w_step):
        wid = width * 10**-3
        thickness = wid/1000
        Re_num = prop[0]*v_max*wid/prop[1]
        Nu_num = Nu_cal(Re_num,prop[3])
        print(Nu_num)
        
        values_3 = {Nu:Nu_num, Re:Re_num, Ts:temp, t:thickness, w:wid}
        eq_Err = eq_Q1.rhs - eq_Q2.rhs
        Error_temp = abs(eq_Err.subs(values_1).subs(values_2).subs(values_3))
        print(Error_temp)
        if Error_temp< Error:
            W = width
            temperature = temp
            Error = Error_temp

print("width = ", end="")
print(W) 
print("Temperature = ", end="")
print(temperature)
print("Error = ", end="")
print(Error)

equation = sp.Eq(eq_Q1.rhs,eq_Q2.rhs)

prop = property(temperature)

values_1 = {l:l_num,s:s_num, Too:Too_num, k:prop[2],w:W*10**-3, t:W*10**-6, Ts:temperature}
eq_1 = equation.subs(values_1)

x = [i for i in range(1,51)]
y = []
print(x)

for v in x:
    Re_num = prop[0]*v*W*10**-3/prop[1]
    Nu_num = Nu_cal(Re_num, prop[3])
    values_2 = {Nu:Nu_num}
    i = sp.solve(eq_1.subs(values_2),I)
    y.append(i[1])

print(y)

plt.plot(x, y)
plt.title('Current vs Velocity')
plt.xlabel('Velocity')
plt.ylabel('Current')
plt.show()
