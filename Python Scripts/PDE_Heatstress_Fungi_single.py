from fipy import *
import numpy as np
from matplotlib import *
import sys
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(threshold=10)
import os
#os.getcwd()
#os.chdir('C:/Users/Felix/Documents/SFB/fungal_interactions')

import pandas as pd

from alternativeViewer import CustomViewer

from alternativeFaceValue import ModifiedCellVariable

######################
#### dimensions ######
######################

nx= 172. # delta_x = delta_y = 0.5 mm; simulate a diameter of 86 mm
ny= nx
dx = 1.
dy = dx
L = nx*dx
 

######################   
#### parameters ######
######################

un= Variable(value=0.97) #substrate uptake of noninsulated biomass
ui = 0.1
c = Variable(value=0.5) #conversion of noninsulated to insulated biomass
#Dn = 10e-7 # diffusion of mobile biomass, m2 s−1
Db = Variable(value=22.6) # diffusion of non-insulated biomass, m2 s−1
beta_1 = 0. # mobilization of biomass
alpha_1 = Variable(value=0.87) # immobilization of biomass

theta = 1. # regulator of switch from net immnobil. to net mobil.
gamma1 = 1. # efficiency of conversion from mobile to immobile biomass (only?)
w = 0.01 #0.01 # replenishment of external substrate
sm = 10. # maximum substrate per cell


# parameters regarding inhibitor production & competition

omega_1 = Variable(value=0.1) # inhibitor production of sp 1
Di1 = Variable(value=10.) # inhibitor diffusion rate of sp 1
chi = 1.0 # metabolic cost of inhibitor production (limiting conversion from n to i)
psi_1 = 0.00000001  # sensitivity of sp 1 to competitor's inhibitor

# change autophagic recycling to degradation term of inhibitor?
# for noe, leave degradation at 0

eta1 = 0.

# parameters regarding heat shock defense

stressbegin = 1.
stresstemp = Variable(value=22.)
endoftreatment= Variable(value=0.)
temperature = Variable(value=22.)
#h_1 = Variable(value=0.) # production rate of hsp for sp1
#h_2 = Variable(value=0.) # production rate of hsp for sp2
hm1 = Variable(value=0.119) # max. production of hsp for sp1
r1 = 0. # decay of hsp
hsplim = Variable(value=0.001)
Dh = 0. # spread of hsp

cutoff = 0.1 # which amount of biomass will be destroyed when heat stress occurs



# when we want to vary parameters, we have tp resets the system after each 
# run, so the loop starts here:
    
from builtins import range

# make dataframe to append following dataframes to
# I use customViewer to make a panda dataframe that I can then save and use in R
        
DF = pd.DataFrame()
    
for i in (10,):
    for j in (0.05,):
        for k in (22,):
        
            # here we reset the parameter to be varied:    
            Db.setValue(i)
            hm1.setValue(j)
            stresstemp.setValue(k)
             
            
            #########################
            ##### produce mesh ######
            #########################
            
            mesh = Grid2D(nx=nx,ny=ny,dx=dx,dy=dy)
                
            x,y = mesh.cellCenters
        
            #####################################
            ##### define state variables ########
            #####################################   
            
        
                    
            # insulated biomass (hyphae with thicker cell walls)
            bi1 = CellVariable(name="bi1",mesh=mesh,hasOld=True,value=0.)
            # noninsulated biomass (tips that grow and take up nutrients)
            bn1 = CellVariable(name="bn1",mesh=mesh,hasOld=True,value=0.)
            # mobile biomass (substrate that has been taken up)
            n1 = CellVariable(name="n1",mesh=mesh,hasOld=True,value=0.)
            # inhibitor
            i1= CellVariable(name="i1",mesh=mesh,hasOld=True,value=0.)   
            # heat shock protein
            hsp1= CellVariable(name="hsp1",mesh=mesh,hasOld=True,value=0.)
            # substrate
            s = CellVariable(name="s",mesh=mesh,hasOld=True,value=0.)
            
            pi1 = CellVariable(name="pi1",mesh=mesh,hasOld=True,value=0.)
          
            # diffusion of mobile biomass n
            Dn1 = ModifiedCellVariable(name="Dn1",mesh=mesh,value=Db*10e-7)
         
            Db1 = ModifiedCellVariable(name="Db1",mesh=mesh,value=Db)
           
        
            
            # bn -> bi conversion parameters as cellvars because inhibitor influences them
            alpha1 = CellVariable(name="alpha1",mesh=mesh,hasOld=True,value=alpha_1)
            beta1 = CellVariable(name="beta1",mesh=mesh,hasOld=True,value=beta_1)
          
            # heat shock protein production rate
            
            h1 = CellVariable(name="h1",mesh=mesh,value=0.)
            ### dummy variable to scale the effect of heat
            z1 = ModifiedCellVariable(name="z1",mesh=mesh,value=1.)
             
            
            
            
            #####################################
            ##### differential equations ########
            #####################################  
            
            ### species 1
            eqbi1 = (TransientTerm(var=bi1)==ImplicitSourceTerm(var=bn1,coeff=c*z1))
            eqbn1 = (TransientTerm(var=bn1)== DiffusionTerm(var=bn1,coeff=Db1.arithmeticFaceValue*z1.onesidedFaceValue)+ImplicitSourceTerm(var=bn1, coeff=(gamma1*z1*alpha1*pi1-z1*beta1*pi1))-ImplicitSourceTerm(var=bn1,coeff=c*z1))
            
            # do i need arithemtic or onesided for Dn1?
            eqn1 = (TransientTerm(var=n1)== DiffusionTerm(var=n1,coeff=Dn1.arithmeticFaceValue*z1.onesidedFaceValue) - ImplicitSourceTerm(var=bn1, coeff=(alpha1*pi1*z1-z1*beta1*pi1))+ImplicitSourceTerm(var=s,coeff=un*z1*bn1)+ImplicitSourceTerm(var=s,coeff=ui*z1*bi1)-ImplicitSourceTerm(var=n1,coeff=omega_1*z1*bn1)-ImplicitSourceTerm(var=n1,coeff=h1*bn1))       
            
            eqi1 = (TransientTerm(var=i1)==ImplicitSourceTerm(var=n1,coeff=omega_1*chi*bn1*z1)+DiffusionTerm(var=i1,coeff=Di1)-ImplicitSourceTerm(var=i1,coeff=eta1))
            eqhsp1 = (TransientTerm(var=hsp1)==ImplicitSourceTerm(var=n1,coeff=h1*bn1)-ImplicitSourceTerm(var=hsp1,coeff=r1)) #+DiffusionTerm(var=hsp1,coeff=Dh)
            
            ## substrate
            eqs = (TransientTerm(var=s)== w*sm-ImplicitSourceTerm(var=s,coeff=w)-ImplicitSourceTerm(var=s,coeff=un*z1*bn1)-ImplicitSourceTerm(var=s,coeff=ui*z1*bi1))
             
            
            
            #####################################
            ######### initial values ############
            #####################################  
            
            bi0 = 0.
            bn0 = 1.
            n0 = 1.
            n_lim = 5. # change? n is max. 6 
            s0 = 10.
            i0= 0.
            r = 6. # radius of initial plug: 3mm = 6 * delta_x
            
            
            #bi1.setValue(bi0,where=((x>(L/2-(r*dx)))&(x<(L/2+(r*dx)))&(y>(L/2-(r*dy)))&(y<(L/2+(r*dy)))))
            #bn1.setValue(bn0,where=((x>(L/2-(r*dx)))&(x<(L/2+(r*dx)))&(y>(L/2-(r*dy)))&(y<(L/2+(r*dy)))))
            #n1.setValue(n0,where=((x>(L/2-(r*dx)))&(x<(L/2+(r*dx)))&(y>(L/2-(r*dy)))&(y<(L/2+(r*dy)))))
            #i1.setValue(i0,where=((x>(L/2-(r*dx)))&(x<(L/2+(r*dx)))&(y>(L/2-(r*dy)))&(y<(L/2+(r*dy)))))
            
            bi1.setValue(bi0,where=(x - L/2.)**2 + (y - L/2)**2 < r**2)
            bn1.setValue(bn0,where=(x - L/2.)**2 + (y - L/2)**2 < r**2)
            n1.setValue(n0,where=(x - L/2.)**2 + (y - L/2)**2 < r**2)
            i1.setValue(i0,where=(x - L/2.)**2 + (y - L/2)**2 < r**2)
            
            s.setValue(s0)
            # heterogeneous substrate:
            #s.setValue(50.,where=(x - L/2.)**2 + (y - L/2)**2 < r**2)
            #s.setValue(50.,where=(x - L/3.)**2 + (y - L/2)**2 < r**2)
            #w = 0.
            #sm = 50.
     
            
            pi1.setValue([x/y if y else 0. for x,y in zip(n1,(bn1+bi1))],where=(bn1+bi1>0)) # /(bn1+bi1)
            Dn1.setValue(Db,where=n1<n_lim)
            
            
            #####################################
            ###### boundary conditions ##########
            #####################################    
            
            valueTopLeft = 0
            valueBottomRight = 0
            X, Y = mesh.faceCenters
            facesTopLeft = ((mesh.facesLeft & (Y > L / 2)) | (mesh.facesTop & (X < L / 2)))
            facesBottomRight = ((mesh.facesRight & (Y < L / 2)) | (mesh.facesBottom & (X > L / 2)))
            
            bi1.constrain(valueTopLeft, facesTopLeft)
            bi1.constrain(valueBottomRight, facesBottomRight)
            bn1.constrain(valueTopLeft, facesTopLeft)
            bn1.constrain(valueBottomRight, facesBottomRight)
            n1.constrain(valueTopLeft, facesTopLeft)
            n1.constrain(valueBottomRight, facesBottomRight)
        
            
            z1.constrain(valueTopLeft, facesTopLeft)
            z1.constrain(valueBottomRight, facesBottomRight)
        
            
            s.constrain(valueTopLeft, facesTopLeft)
            s.constrain(valueBottomRight, facesBottomRight)
        
        
        
        
            
            #-------------------------------------
            
            ########################
            ###### simulation ######
            ########################
            
                
            eq = eqbi1 & eqbn1 & eqn1 & eqs & eqi1 & eqhsp1
            
            vi = Viewer(vars=(bi1+bn1+n1,i1))
        
        
            
            for loopvar in range(0,101):
        
                t = loopvar / 10.
                
                # 1t = 1d -> 1h = 0.04t
                
                print(t)                       
                # update pi
                
                ##########################
                ###### heat effects ######
                ##########################
                
                # turn heat stress on/off:
                if t == stressbegin:
                    temperature.setValue(stresstemp)
                    if temperature >= 45.:
                        z1.setValue(0.,where=(n1+bi1+bn1)>=cutoff)
                        bi1.setValue(0.,where=z1==1.)
                        bn1.setValue(0.,where=z1==1.)
                        n1.setValue(0.,where=z1==1.)
        
                        h1.setValue((hm1*((temperature-22)/(45-22))),where=n1+bi1+bn1>=cutoff)
                        
                if t == stressbegin+0.1: # ca. 2h later
                    temperature.setValue(22.)
                    endoftreatment.setValue(t)
                
                h1.setValue(0.,where=((hsp1>hsplim*(bi1+bn1+n1))))
                z1.setValue(1.,where=((hsp1>=hsplim*(bi1+bn1+n1))))
        
        
                #################################
                #### no competition effects #####
                #################################
                
                pi1.setValue([x/y if y else 0. for x,y in zip(n1,(bn1+bi1))],where=(bn1+bi1>0)) # /(bn1+bi1)
                pi1.setValue(0.,where=(bn1+bi1<=0))
                pi1.setValue(0.,where=(n1<=0))
        
                Dn1.setValue(Db1,where=n1<n_lim)
                Dn1.setValue(Db1*10e-7,where=n1>n_lim)
                
                 
                 # update all PDEs
                bi1.updateOld()
                bn1.updateOld()
                n1.updateOld()
                i1.updateOld()
                hsp1.updateOld()
                s.updateOld()
                
                if (t == 10):
                     x = CustomViewer(vars=(bn1,bi1,n1,pi1,Dn1,z1,hsp1,i1,h1,s)).makedf(time=t)
                     x.insert(2,"Db",i)
                     x.insert(2,"hm",j)
                     x.insert(2,"temperature",k)
                     DF = DF.append(x)
                     #TSVViewer(vars=(bn1,bi1,n1,pi1,Dn1,z1,hsp1,bn2,bi2,n2,pi2,Dn2,z2,hsp2,s)).plot(filename="dual_hsp_t{0}.txt".format(t))
                     vi = Viewer(vars=(n1+bi1+bn1,i1),FIPY_VIEWER="matplotlib")
                 
                eq.solve(dt=0.1)
                
        
            # save final dataframe 
    
DF.to_csv("testrun.txt",index=False)
