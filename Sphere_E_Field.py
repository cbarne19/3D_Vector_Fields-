#Sphere 
import numpy as np 
import pandas as pd 
import plotly.graph_objects as go

x = np.linspace(-20,20,num=1000)
y = np.linspace(-20,20,num=1000)

Q = 5 

def Hollow_Sphere_maker(R,res,charge):
    data = [] 
    xsize = 2*R 
    ysize = 2*R  
    zsize = 2*R
    for i in range(res):
        for j in range(res):
            for k in range(res):
                x = (i*(xsize/res)-xsize/2)
                y = (j*(xsize/res)-ysize/2)
                z = (k*(xsize/res)-zsize/2)
                if (x)**2 + y**2 +z**2 <= R**2 and (x)**2 + y**2 +z**2 >= (R-0.01)**2:
                    data.append([x,y,z,charge])
                else:
                    continue  
    return data


def pythag(x1,x2,y1,y2):
    return ((abs(x1-x2)**2)+(abs(y1-y2)**2))**0.5 

def pythag3d(x1,x2,y1,y2,z1,z2): 
    return ((abs(x1-x2)**2)+(abs(y1-y2)**2)+(abs(z1-z2))**2)**0.5 

def E_field_calc(x_c,y_c,z_c,Q,qs,n):
    E_0 = 8.85e-12
    k = 1/(4*np.pi*E_0)
    xsize = int(5*max(x_c)) 
    ysize = int(5*max(y_c)) 
    zsize = int(5*max(z_c))
    Tot_E_Field = [] 
    data = [] 
    for i in range(n+1):
        for j in range(n+1): 
            for m in range(n+1):
                x = i*(xsize/(n+1))-(xsize/2)
                y = j*(ysize/(n+1))-(ysize/2)
                z = m*(zsize/(n+1))-(zsize/2)
                for l in range(len(x_c)):
                    summ = (eval(qs[l])/len(x_c))/((pythag3d(x,x_c[l],y,y_c[l],z,z_c[l]))**2) 
                    unit_length = pythag3d(x,x_c[l],y,y_c[l],z,z_c[l])
                    unit_vector = (np.array([x,y,z])-np.array([x_c[l],y_c[l],z_c[l]]))/unit_length 
                    E_field_point = k*summ*unit_vector
                    Tot_E_Field.append(E_field_point)
                Summed_E = sum(Tot_E_Field) 
                data.append([x,y,z,Summed_E[0],Summed_E[1],Summed_E[2]])
                Tot_E_Field = []             
    return data

R = 5 
res = 100 

data2 = Hollow_Sphere_maker(R,res,'Q')

df = pd.DataFrame(data2, columns = ['x', 'y','z','charge'])

xs = np.ndarray.tolist(np.array(df['x']))
ys = np.ndarray.tolist(np.array(df['y']))
zs = np.ndarray.tolist(np.array(df['z']))
qs = np.ndarray.tolist(np.array(df['charge']))

data3 = E_field_calc(xs,ys,zs,Q,qs,10)

df2 = pd.DataFrame(data3, columns = ['x', 'y','z','u','v','w'])

fig = go.Figure(data = go.Cone(
    x=df2['x'],
    y=df2['y'],
    z=df2['z'],
    u=df2['u'], 
    v=df2['v'],  
    w=df2['w'], 
    colorscale='Blues',
    sizemode="scaled",
    sizeref=1))

fig.add_trace(go.Scatter3d(
    x=df["x"], y=df["y"], z=df["z"],
    mode = 'markers'))

fig.update_layout(scene=dict(aspectratio=dict(x=1, y=1, z=1),
                         camera_eye=dict(x=1, y=1, z=1)))

fig.write_html("Sphere_field.html") 