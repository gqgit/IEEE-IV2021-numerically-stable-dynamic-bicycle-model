"""
Numerically Stable Dynamic Bicycle Model for Discrete-time Control (IEEE IV2021 workshop paper)
this code is the demonstration of numerical stability of our method versus forward-Euler discretization.
@Author: Qiang Ge
@Date: 2021.06.03
@email: gq17@mails.tsinghua.edu.cn
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

T  = 0.01
T_benchmark = 0.001
STEPS = int(4/T)

speed = 8

Rw = 0.325       #rolling radius of drive wheel, m
m = 1412.0       #mass of vehicle, kg
Iz = 1536.7      #yaw inertia of vehicle, kg*m^2
a = 1.06         #distance between C.G. and front axle, m
b = 1.85         #distance between C.G. and rear axle, m
kf = -128915.5   #equivalent sideslip force stiffness of front axle, N/rad
kr = -85943.6    #equivalent sideslip force stiffness of rear axle, N/rad


FONT_SIZE = 18
FONT_SIZE_ = 14
family = 'Times New Roman'

def dynamic_linear_backward(x0,u0,T):
    '''
    dynamic 3DoF bicycle model, analytical, defferentiable.
    linear tire sideslip force, with stiffness estimation: kf and kr.
    the discretization is a tailored method inspired by (but not exactly be) backward Euler method.
    '''
    
    x1 = np.zeros(len(x0))
    
    x1[0] = x0[0] + T*(x0[3]*np.cos(x0[2])-x0[4]*np.sin(x0[2]))
    x1[1] = x0[1] + T*(x0[4]*np.cos(x0[2])+x0[3]*np.sin(x0[2]))
    x1[2] = x0[2] + T* x0[5]
    x1[3] = x0[3] + T* u0[0]
    x1[4] = (-(a*kf-b*kr)*x0[5]+kf*u0[1]*x0[3]+m*x0[5]*x0[3]*x0[3]-m*x0[3]*x0[4]/T)/(kf+kr-m*x0[3]/T)
    x1[5] = (-Iz*x0[5]*x0[3]/T-(a*kf-b*kr)*x0[4]+a*kf*u0[1]*x0[3])/((a*a*kf+b*b*kr)-Iz*x0[3]/T)
    
    return x1


def dynamic_linear_forward(x0,u0,T):
    '''
    dynamic 3DoF bicycle model, analytical, defferentiable.
    linear tire sideslip force, with stiffness estimation: kf and kr.
    the discretization is via standard forward Euler method.
    drawback: this model could be numerically instable (e.g. when T = 0.1 sec)
              and has a singularity issue when speed approaches zero.
    '''
    
    x1 = np.zeros(len(x0))

    x1[0] = x0[0] + T*(x0[3]*np.cos(x0[2])-x0[4]*np.sin(x0[2]))
    x1[1] = x0[1] + T*(x0[4]*np.cos(x0[2])+x0[3]*np.sin(x0[2]))
    x1[2] = x0[2] + T* x0[5]
    x1[3] = x0[3] + T*(u0[0] + x0[4]*x0[5] - kf*(x0[4]+a*x0[5]-u0[1]*x0[3])*np.sin(u0[1])/m/x0[3])
    x1[4] = x0[4] + T*(kf*(x0[4]+a*x0[5]-u0[1]*x0[3])+kr*(x0[4]-b*x0[5])-m*x0[3]*x0[3]*x0[5])/m/x0[3]
    x1[5] = x0[5] + T*(a*kf*(x0[4]+a*x0[5]-u0[1]*x0[3])-b*kr*(x0[4]-b*x0[5]))/Iz/x0[3]
    
    return x1


def trajectory(x0,ur,T,N,type):
    #used for taylor expansion and linearization!
    traj = []
    traj.append(x0)

    if type == 'linear_forward':
        for i in range(N):
           result = dynamic_linear_forward(traj[i],ur[i],T)
           traj.append(result)
    else:
        for i in range(N):
           result = dynamic_linear_backward(traj[i],ur[i],T)
           traj.append(result)

    return traj

if __name__ == "__main__":
    
    x0 = [0.,0.,0.,speed,0.,0.]

####step input
    '''input: desired acceleration [m/s^2] and front wheel steer angle [rad]
              desired acceleration = driving torque/m/Rw
    '''
    ur_back = [[0.,0.2674]]*STEPS
    ur_forw = [[0.,0.2674]]*STEPS
    
####sine input
#    ur_back = []
#    ur_forw = []
#    for i in range(STEPS):
#        ur_back.append([0.,0.3 * np.sin(np.pi*i*T)])
#        ur_forw.append([0.,0.3 * np.sin(np.pi*i*T)])
    

###simulink step response###
    benchmark = loadmat("data/simulink_step_u_equals_"+str(speed)+".mat")
#    benchmark = loadmat("data/simulink_sine_u_equals_"+str(speed)+".mat")
    w_benchmark = benchmark["omega"].tolist()
    v_benchmark = benchmark["v"].tolist()
    t_benchmark = [i*T_benchmark for i in list(range(len(w_benchmark)))]
    benchmark_label = 'Simulink'
    y_lim_v = [-0.05,1.25]
    y_lim_w = [-0.05,1.25]


    tra_f = trajectory(x0,ur_forw,T,STEPS,'linear_forward')
    tra_b = trajectory(x0,ur_back,T,STEPS,'linear_backward')

    x_f = [];   x_b = []
    y_f = [];   y_b = []
    yaw_f = []; yaw_b = []
    u_f = [];   u_b = []
    v_f = [];   v_b = []
    w_f = [];   w_b = []

    for i in range(STEPS+1):
        x_f.append(tra_f[i][0])
        y_f.append(tra_f[i][1])
        yaw_f.append(tra_f[i][2])
        u_f.append(tra_f[i][3])
        v_f.append(tra_f[i][4])
        w_f.append(tra_f[i][5])
        
        x_b.append(tra_b[i][0])
        y_b.append(tra_b[i][1])
        yaw_b.append(tra_b[i][2])
        u_b.append(tra_b[i][3])
        v_b.append(tra_b[i][4])
        w_b.append(tra_b[i][5])

    
    
    F2 = plt.figure(figsize=(8,4),dpi=100)
    
    ax9 = F2.add_subplot(1,2,1)
    ax9.plot([i * T for i in range(STEPS+1)],v_f,'-r',label='forward Euler')
    ax9.plot([i * T for i in range(STEPS+1)],v_b,'-b',label='our method')
    ax9.plot(t_benchmark, v_benchmark,'-g',label=benchmark_label)
    ax9.set_xlabel('time (s)', fontdict={'family':family,'size': FONT_SIZE})
    ax9.set_ylabel('lateral velocity (m/s)', fontdict={'family':family,'size': FONT_SIZE} )
    ax9.set_ylim(y_lim_v)
    plt.xticks(fontproperties = family, size = FONT_SIZE_)
    plt.yticks([0.0,0.4,0.8,1.2],fontproperties = family, size = FONT_SIZE_)
    ax9.legend(loc='center right',prop={'family': family,'size': FONT_SIZE_})
    
    ax10 = F2.add_subplot(1,2,2)
    ax10.plot([i * T for i in range(STEPS+1)],w_f,'-r',label='forward Euler')
    ax10.plot([i * T for i in range(STEPS+1)],w_b,'-b',label='our method')
    ax10.plot(t_benchmark, w_benchmark,'-g',label=benchmark_label)
    ax10.set_xlabel('time (s)', fontdict={'family':family,'size': FONT_SIZE})
    ax10.set_ylabel('yaw rate (rad/s)', fontdict={'family':family,'size': FONT_SIZE});
    ax10.set_ylim(y_lim_w)
    plt.xticks(fontproperties = family, size = FONT_SIZE_)
    plt.yticks([0.0,0.4,0.8,1.2],fontproperties = family, size = FONT_SIZE_)
    legend = ax10.legend(loc='upper right',prop={'family': family,'size': FONT_SIZE_})

    
    plt.tight_layout()
    
    plt.show()
    
    
    
    
    
    
    
