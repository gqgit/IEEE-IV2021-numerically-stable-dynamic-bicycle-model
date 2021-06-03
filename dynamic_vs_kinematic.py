"""
Numerically Stable Dynamic Bicycle Model for Discrete-time Control (IEEE IV2021 workshop paper)
this code is the demonstration of simulation accuracy of our method versus kinematic model.
@Author: Qiang Ge
@Date: 2021.06.03
@email: gq17@mails.tsinghua.edu.cn
"""

import numpy as np
import matplotlib.pyplot as plt
import csv

T  = 0.1
STEPS = int(4/T)+1

speed = 8

Rw = 0.325       #rolling radius of drive wheel, m
m = 1412.0       #mass of vehicle, kg
Iz = 1536.7      #yaw inertia of vehicle, kg*m^2
a = 1.06         #distance between C.G. and front axle, m
b = 1.85         #distance between C.G. and rear axle, m
kf = -128915.5   #equivalent sideslip force stiffness of front axle, N/rad
kr = -85943.6    #equivalent sideslip force stiffness of rear axle, N/rad

family = 'Times New Roman'
FONT_SIZE  = 18
FONT_SIZE_ = 14

def dynamic_linear_backward(x0,u0,T):
    '''
    dynamic 3DoF bicycle model, analytical, defferentiable.
    linear tire sideslip force, with stiffness estimation: kf and kr.
    the discretization is a tailored method inspired (but not exactly) backward Euler method.
    '''
    
    x1 = np.zeros(len(x0))
    
    x1[0] = x0[0] + T*(x0[3]*np.cos(x0[2])-x0[4]*np.sin(x0[2]))
    x1[1] = x0[1] + T*(x0[4]*np.cos(x0[2])+x0[3]*np.sin(x0[2]))
    x1[2] = x0[2] + T* x0[5]
    x1[3] = x0[3] + T* u0[0]
    x1[4] = (-(a*kf-b*kr)*x0[5]+kf*u0[1]*x0[3]+m*x0[5]*x0[3]*x0[3]-m*x0[3]*x0[4]/T)/(kf+kr-m*x0[3]/T)
    x1[5] = (-Iz*x0[5]*x0[3]/T-(a*kf-b*kr)*x0[4]+a*kf*u0[1]*x0[3])/((a*a*kf+b*b*kr)-Iz*x0[3]/T)
    
    return x1

def kinematic_forward(x0,u0,T):
    '''
    kinematic 2DoF bicycle model, analytical, defferentiable.
    the discretization is via standard forward Euler method.
    '''

    x1 = np.zeros(len(x0))

    x1[0] = x0[0] + T*(x0[3]*np.cos(x0[2]) - b*np.sin(x0[2])*x0[3]*np.tan(u0[1])/(a+b))
    x1[1] = x0[1] + T*(x0[3]*np.sin(x0[2]) + b*np.cos(x0[2])*x0[3]*np.tan(u0[1])/(a+b))
    x1[2] = x0[2] + T*(x0[3]*np.tan(u0[1])/(a+b))
    x1[3] = x0[3] + T* u0[0]
    x1[4] = 0.
    x1[5] = 0.
    return x1


def trajectory(x0,ur,T,N,type):
    traj = []
    traj.append(x0)

    if type == 'kinematic':
        for i in range(N):
            result = kinematic_forward(traj[i],ur[i],T)
            traj.append(result)
    elif type == 'linear_backward':
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
    ur = [[0.,0.2674]]*STEPS #step input


    x_carsim_str = []
    y_carsim_str = []
    v_carsim_str = []
    w_carsim_str = []
    t_carsim_str = []
    string_speed = str(speed)
    string_speed_no_dot = ''
    for i in string_speed:
        if i != '.':
            string_speed_no_dot = string_speed_no_dot + i
    csv_reader = csv.reader(open("data/carsim_u_equals_"+string_speed_no_dot+".csv"))

    kk = 0
    for i in csv_reader:
        if kk > 0:
            t_carsim_str.append(i[0])   #time, 1 column, s
            x_carsim_str.append(i[551]) #x, 552 column, m
            y_carsim_str.append(i[571]) #vx, 572 column, kph
            v_carsim_str.append(i[525]) #vy, 526 column, kph
            w_carsim_str.append(i[34])  #omega, 35 column, deg/s
        kk = kk + 1
    t_carsim = list(map(eval,t_carsim_str))
    x_carsim = list(map(eval,x_carsim_str))
    y_carsim = list(map(eval,y_carsim_str))
    v_carsim = list(map(eval,v_carsim_str))
    w_carsim = list(map(eval,w_carsim_str))

    tra_f = trajectory(x0,ur,T,STEPS,'kinematic')
    tra_b = trajectory(x0,ur,T,STEPS,'linear_backward')

    x_f = [];   x_b = []
    y_f = [];   y_b = []
    yaw_f = []; yaw_b = []
    u_f = [];   u_b = []
    w_b = []
    v_b = []

    for i in range(STEPS+1):
        x_f.append(tra_f[i][0])
        y_f.append(tra_f[i][1])
        yaw_f.append(tra_f[i][2])
        u_f.append(tra_f[i][3])
        
        x_b.append(tra_b[i][0])
        y_b.append(tra_b[i][1])
        yaw_b.append(tra_b[i][2])
        u_b.append(tra_b[i][3])
        v_b.append(tra_f[i][4])
        w_b.append(tra_f[i][5])

        
    time_prev = []
    err_forw = []; err_forw_rms = []
    err_back = []; err_back_rms = []
    for i in range(STEPS):
        for j in range(len(t_carsim)):
            if abs(i*T - t_carsim[j]) < 0.0001:#numerical accuracy issue of data, no special meaning
                time_prev.append(i*T)
                err_forw.append( ((x_carsim[j]+1.058-x_f[i])**2+(y_carsim[j]-y_f[i])**2) )
                err_back.append( ((x_carsim[j]+1.058-x_b[i])**2+(y_carsim[j]-y_b[i])**2) )
                err_forw_rms.append( (sum(err_forw)/len(err_forw))**0.5 )
                err_back_rms.append( (sum(err_back)/len(err_back))**0.5 )

    
    F1 = plt.figure(figsize=(4.5,4.5),dpi=100)
    
    ax1 = F1.add_subplot(1,1,1)
    ax1.plot(x_f,y_f,'-r',label='kinematic (forward Euler)')
    ax1.plot(x_b,y_b,'-b',label='dynamic (our method)')
    ax1.plot([i+1.058 for i in x_carsim],y_carsim,'-g',label='CarSim')#+1.058 because CarSim reference point is at front axle instead of Center of Gravity (C.G.)
    ax1.set_xlabel('x location (m)', fontdict={'family':family,'size': FONT_SIZE} )
    ax1.set_ylabel('y location (m)', fontdict={'family':family,'size': FONT_SIZE} )
    plt.xticks(fontproperties = 'Times New Roman', size = FONT_SIZE_)
    plt.yticks(fontproperties = 'Times New Roman', size = FONT_SIZE_)
    ax1.legend(loc='center left',prop={'family': family,'size': FONT_SIZE_}, framealpha = 0.5)
    ax1.axis('equal')
    ax1.axis([-10,13,-1,23])
    plt.tight_layout()

    F2 = plt.figure(figsize=(4.5,4.5),dpi=100)
    ax11 = F2.add_subplot(1,1,1)
    ax11.plot(time_prev,err_forw_rms,'-r',label='kinematic (forward Euler)')
    ax11.plot(time_prev,err_back_rms,'-b',label='dynamic (our method)')
    ax11.set_xlabel('prediction time (s)', fontdict={'family':family,'size': FONT_SIZE} )
    plt.xticks(fontproperties = 'Times New Roman', size = FONT_SIZE_)
    plt.yticks([0,1,2,3,4],fontproperties = 'Times New Roman', size = FONT_SIZE_)
    ax11.set_ylabel('predictive rms error (m)', fontdict={'family':family,'size': FONT_SIZE} )
    ax11.legend(loc='upper left',prop={'family': family,'size': FONT_SIZE_}, framealpha = 0.5)
    
    plt.tight_layout()
    
    print('spped:       ',speed)
    print('err_dyna:    ',err_back_rms[-1])
    print('err_kine:    ',err_forw_rms[-1])
    print('improved by: ',(err_forw_rms[-1] - err_back_rms[-1])/err_forw_rms[-1])
    
    plt.show()
    
    
        
        
    
    
    
    
    
    
    