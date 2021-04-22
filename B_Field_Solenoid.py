import numpy as np 
from numpy import sin,pi,cos
import pandas as pd 
import plotly.graph_objects as go

q = 5 

def pythag(x1,x2,y1,y2):
    return ((abs(x1-x2)**2)+(abs(y1-y2)**2))**0.5 

def pythag3d(x1,x2,y1,y2,z1,z2): 
    return ((abs(x1-x2)**2)+(abs(y1-y2)**2)+(abs(z1-z2))**2)**0.5 
             
def B_field_calc(x_c,y_c,z_c,I,n):
    mu_0 = 4*pi*10**(-7)
    k = mu_0/(4*np.pi)
    xsize = 5*max(x_c) 
    ysize = 5*max(y_c) 
    zsize = 5*max(z_c)
    Tot_B_Field = [] 
    data = [] 
    for i in range(n+1):
        for j in range(n+1): 
            for m in range(n+1):
                x = i*(xsize/(n+1))-(xsize/2)
                y = j*(ysize/(n+1))-(ysize/2)
                z = m*(zsize/(n+1))-(zsize/2)
                for l in range(len(x_c)-1):                        
                    summ = (I*pythag3d(x_c[l+1],x_c[l],y_c[l+1],y_c[l],z_c[l+1],z_c[l]))/((pythag3d(x,x_c[l],y,y_c[l],z,z_c[l]))**2) 
                    unit_length = pythag3d(x,x_c[l],y,y_c[l],z,z_c[l])
                    unit_len_dl = pythag3d(x_c[l+1],x_c[l],y_c[l+1],y_c[l],z_c[l+1],z_c[l])
                    unit_vector = (np.array([x,y,z])-np.array([x_c[l],y_c[l],z_c[l]]))/unit_length 
                    dl_unit = (np.array([x_c[l],y_c[l],z_c[l]])-np.array([x_c[l+1],y_c[l+1],z_c[l+1]]))/unit_len_dl
                    B_field_point = k*summ*np.cross(unit_vector,dl_unit)
                    Tot_B_Field.append(B_field_point)
                Summed_B = sum(Tot_B_Field) 
  
                data.append([x,y,z,Summed_B[0],Summed_B[1],Summed_B[2],'']) 
                Tot_B_Field = []      
    return data

def solenoid_maker(L,N,r):
    x = np.linspace(-L/2,L/2,num=5000) 
    data = []
    for i in range(len(x)):
        data.append([x[i],r*sin(((2*pi*N)/L)*x[i]),
                     r*cos(((2*pi*N)/L)*x[i]),'Solenoid'])
    return data 

data2 = solenoid_maker(15e-2,30,4.2e-2)


df = pd.DataFrame(data2, columns = ['x', 'y','z','Name'])

xs = np.ndarray.tolist(np.array(df['x']))
ys = np.ndarray.tolist(np.array(df['y']))
zs = np.ndarray.tolist(np.array(df['z']))

I = 16

data3 = B_field_calc(xs,ys,zs,I,10)

df2 = pd.DataFrame(data3, columns = ['x', 'y','z','u','v','w','Angle'])


fig = go.Figure(data = go.Cone(
    x=df2['x'],
    y=df2['y'],
    z=df2['z'],
    u=df2['u'], 
    v=df2['v'],  
    w=df2['w'], 
    colorscale='Blues',
    sizemode="scaled",
    sizeref=1.5))

fig.add_trace(go.Scatter3d(
    x=df["x"], y=df["y"], z=df["z"],
    marker=dict(
        size=0,
        color='#ff0000',
        colorscale='Viridis',
    ),
    line=dict(
        color='#ff0000',
        width=10
    )
))
     
fig.update_layout(scene=dict(aspectratio=dict(x=1, y=1, z=1),
                         camera_eye=dict(x=1, y=1, z=1)))

fig.update_layout(title_text="Magnetic vector Field for 15cm 30 Turn solenoid with current of 16A")

fig.update_layout(scene = dict(
                    xaxis_title='x (m)',
                    yaxis_title='y (m)',
                    zaxis_title='z (m)') ,template = 'plotly_dark'
                   
                   )
fig.write_html("Solenoid_Test_Biot_W_FIELD_4.html") 