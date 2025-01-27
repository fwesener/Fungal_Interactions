# -*- coding: utf-8 -*-

from fipy import *
import numpy as np
from matplotlib import *
import sys
#np.set_printoptions(threshold=sys.maxsize) # this command changes how many decimals of a number are shown in the console
#np.set_printoptions(threshold=10)
import os
os.getcwd()
#os.chdir('C:/Users/Felix/Documents/SFB/fungal_interactions')
import pandas as pd
import numpy as np

# the following packages I created to avoid some issues with diffusion at the colony periphery
from alternativeViewer import CustomViewer
from alternativeFaceValue import ModifiedCellVariable

######################
#### dimensions ######
######################

dx = 1.# delta_x = delta_y = 1; one grid cell is 0.5 mm long;
dy = dx
nx= 172. # simulate a diameter of 86 mm = 172 grid cells
ny= nx
L = nx*dx # length of one side of the grid
 
######################   
#### parameters ######
######################

# some paramters differ over space and are thus declared as cell variables further down; 
# these cell variables have names without _, while their numeric values are values marked with _
# if values differ for species 1 and 2, they are marked with numbers

un= 0.97 # substrate uptake of noninsulated biomass
ur = 0.1 # substrate uptake of rigid, insulated biomass
c = Variable(value=0.5) #conversion of noninsulated to insulated biomass
#Dn = 10e-7 # diffusion of mobile biomass, m2 s−1
Dn_1 = Variable(value=22.6) # diffusion of non-insulated biomass of species 1, m2 s−1
Dn_2 = Variable(value=1.005) # diffusion of non-insulated biomass of species 2, m2 s−1
beta_1 = 0. # mobilization of biomass
alpha_1 = 0.87 # immobilization of biomass
alpha_2 = 0.87
beta_2 = 0.

gamma1 = 1. # efficiency of conversion from mobile to immobile biomass 
w = 0.01  # replenishment of external substrate
sm = 10. # maximum substrate per cell

##############################################
##### parameters regarding inhibitor production & competition

omega_1 = Variable(value=0.1) # inhibitor production 
omega_2 = Variable(value=0.1) 
Di1 = Variable(value=0.) # inhibitor diffusion rate
Di2 = Variable(value=0.) 
chi = 1. # metabolic cost of inhibitor production (limiting conversion from n to i)
psi_1 = Variable(value=0.00000001)  # resistance of sp 1 to competitor's inhibitor
psi_2 = Variable(value=0.00000001)  # resistance of sp 2 to competitor's inhibitor

# degradation term of inhibitor?
# for now, leave degradation at 0

eta1 = 0.
eta2 = 0.

##############################################
##### parameters regarding heat shock defense

temperature = Variable(value=22.) # control temperature

stressbegin = 1.
stresstemp = Variable(value=22.) # 22 for control conditions, 45 for stress treatment


hm1 = Variable(value=0.119) # max. production of hsp for species 1
hm2 = Variable(value=0.047) # max. production of hsp for species 2
r1 = 0. # decay of hsp
r2 = 0.
hsplim = Variable(value=0.001)

cutoff = 0.1 # which amount of biomass will be destroyed when heat stress occurs


# when we want to vary parameters, we have tp resets the system after each 
# run, so the loop starts here:
    
from builtins import range

# make dataframe to append following dataframes to
# I use customViewer to make a panda dataframe that I can then save and use in R
        
DF = pd.DataFrame()

# new: vary (focal) species 1 discreetly (three values of growth and lag phase)
# vary species 2 continuously
# leave omega out but vary Di yes/no
for i in (5,17.5):
    for j in range(1,22,2):
        for k in (0.05,0.02): # 0.05,0.02
            for l in np.arange(0.01,0.2,0.02): #0.000001,
                for m in (0,10):
                    for q in (0,10):
                        for te in (22.,45.):
                                # here we reset the parameter to be varied: 
                                Dn_1.setValue(i)
                                Dn_2.setValue(j)
                                
                                hm1.setValue(k)
                                hm2.setValue(l)
                                Di1.setValue(m)
                                Di2.setValue(q)
                                
                                stresstemp.setValue(te)
                                # set psi to assess effect of resistance values:
                                #psi_1.setValue(m)
                                #psi_2.setValue(m)
                                
                            
                            
                                #########################
                                ##### produce mesh ######
                                #########################
                            
                                mesh = Grid2D(nx=nx,ny=ny,dx=dx,dy=dy)
                                
                                x,y = mesh.cellCenters
                                
                                #####################################
                                ##### define state variables ########
                                #####################################    
                                    
                                # insulated biomass (rigid hyphae with thicker cell walls)
                                br1 = CellVariable(name="br1",mesh=mesh,hasOld=True,value=0.)
                                br2 = CellVariable(name="br2",mesh=mesh,hasOld=True,value=0.)
                                
                                # noninsulated biomass (tips that grow and take up nutrients and are not yet rigidified)
                                bn1 = CellVariable(name="bn1",mesh=mesh,hasOld=True,value=0.)
                                bn2 = CellVariable(name="bn2",mesh=mesh,hasOld=True,value=0.)
                                
                                # mobile biomass (substrate that has been taken up)
                                bm1 = CellVariable(name="bm1",mesh=mesh,hasOld=True,value=0.)
                                bm2 = CellVariable(name="bm2",mesh=mesh,hasOld=True,value=0.)
                                
                                # inhibitor
                                i1= CellVariable(name="i1",mesh=mesh,hasOld=True,value=0.)
                                i2= CellVariable(name="i2",mesh=mesh,hasOld=True,value=0.)
                                
                                # heat shock protein
                                hsp1= CellVariable(name="hsp1",mesh=mesh,hasOld=True,value=0.)
                                hsp2= CellVariable(name="hsp2",mesh=mesh,hasOld=True,value=0.)
                                
                                
                                pi1 = CellVariable(name="pi1",mesh=mesh,hasOld=True,value=0.)
                                pi2 = CellVariable(name="pi2",mesh=mesh,hasOld=True,value=0.)
                                
                                # diffusion of mobile biomass bm
                                Dm1 = ModifiedCellVariable(name="Dm1",mesh=mesh,value=Dn_1*10e-7)
                                Dm2 = ModifiedCellVariable(name="Dm2",mesh=mesh,value=Dn_2*10e-7)
                                
                                # diffusion of hyphae bn
                                Dn1 = ModifiedCellVariable(name="Dn1",mesh=mesh,value=Dn_1)
                                Dn2 = ModifiedCellVariable(name="Dn2",mesh=mesh,value=Dn_2)
                                
                                # substrate
                                s = CellVariable(name="s",mesh=mesh,hasOld=True,value=0.)
                                
                                # bn -> bi conversion parameters as cellvars because inhibitor influences them
                                alpha1 = CellVariable(name="alpha1",mesh=mesh,hasOld=True,value=alpha_1)
                                beta1 = CellVariable(name="beta1",mesh=mesh,hasOld=True,value=beta_1)
                                alpha2 = CellVariable(name="alpha2",mesh=mesh,hasOld=True,value=alpha_2)
                                beta2 = CellVariable(name="beta2",mesh=mesh,hasOld=True,value=beta_2)
                                
                                # heat shock protein production rate
                                
                                h1 = CellVariable(name="h1",mesh=mesh,value=0.)
                                h2 = CellVariable(name="h2",mesh=mesh,value=0.)
                                
                                ### dummy variable to scale the effect of heat
                                z1 = ModifiedCellVariable(name="z1",mesh=mesh,value=1.)
                                z2 = ModifiedCellVariable(name="z2",mesh=mesh,value=1.)
                                
                                
                                # Db,Dn onesided or arithmetic?
                                
                                #####################################
                                ##### differential equations ########
                                #####################################    
                                
                                ### species 1
                                
                                # rigid biomass br (old hyphae)
                                eqbr1 = (TransientTerm(var=br1)==ImplicitSourceTerm(var=br1,coeff=c*z1))
                                # noninsukated biomass bn (tips & young hyphae)
                                eqbn1 = (TransientTerm(var=bn1)==DiffusionTerm(var=bn1,coeff=Dn1.arithmeticFaceValue*z1.onesidedFaceValue)+ImplicitSourceTerm(var=bn1, coeff=(gamma1*z1*alpha1*pi1-z1*beta1*pi1))-ImplicitSourceTerm(var=bn1,coeff=c*z1))
                                # mobile biomass bm (internal substrate)
                                eqbm1 = (TransientTerm(var=bm1)== DiffusionTerm(var=bm1,coeff=Dm1.arithmeticFaceValue*z1.onesidedFaceValue) - ImplicitSourceTerm(var=bn1, coeff=(alpha1*pi1*z1-z1*beta1*pi1))+ImplicitSourceTerm(var=s,coeff=un*z1*bn1)+ImplicitSourceTerm(var=s,coeff=ur*z1*br1)-ImplicitSourceTerm(var=bm1,coeff=omega_1*z1*bn1)-ImplicitSourceTerm(var=bm1,coeff=h1*bn1))       
                                # inhibitor i
                                eqi1 = (TransientTerm(var=i1)==ImplicitSourceTerm(var=bm1,coeff=omega_1*chi*bn1*z1)+DiffusionTerm(var=i1,coeff=Di1)-ImplicitSourceTerm(var=i1,coeff=eta1))
                                #heat shock protein hsp
                                eqhsp1 = (TransientTerm(var=hsp1)==ImplicitSourceTerm(var=bm1,coeff=h1*bn1)-ImplicitSourceTerm(var=hsp1,coeff=r1)) #+DiffusionTerm(var=hsp1,coeff=Dh)
                                
                                ### species 2
                                eqbr2 = (TransientTerm(var=br2)==ImplicitSourceTerm(var=bn2,coeff=c*z2))
                                eqbn2 = (TransientTerm(var=bn2)== DiffusionTerm(var=bn2,coeff=Dn2.arithmeticFaceValue*z2.onesidedFaceValue)+ImplicitSourceTerm(var=bn2, coeff=(gamma1*z2*alpha2*pi2-z2*beta2*pi2))-ImplicitSourceTerm(var=bn2,coeff=c*z2))
                                eqbm2 = (TransientTerm(var=bm2)== DiffusionTerm(var=bm2,coeff=Dm2.arithmeticFaceValue*z2.onesidedFaceValue) - ImplicitSourceTerm(var=bn2, coeff=(alpha2*pi2*z2-z2*beta2*pi2))+ImplicitSourceTerm(var=s,coeff=un*z2*bn2)+ImplicitSourceTerm(var=s,coeff=ur*z2*br2)-ImplicitSourceTerm(var=bm2,coeff=omega_2*z2*bn2)-ImplicitSourceTerm(var=bm2,coeff=h2*bn2))       
                                eqi2 = (TransientTerm(var=i2)==ImplicitSourceTerm(var=bm2,coeff=omega_2*chi*bn2*z2)+DiffusionTerm(var=i2,coeff=Di2)-ImplicitSourceTerm(var=i2,coeff=eta2))
                                eqhsp2 = (TransientTerm(var=hsp2)==ImplicitSourceTerm(var=bm2,coeff=h2*bn2)-ImplicitSourceTerm(var=hsp2,coeff=r2)) #+DiffusionTerm(var=hsp1,coeff=Dh)
                                
                                ## substrate is depleted by species 1 and 2, and can be replenished via w
                                eqs = (TransientTerm(var=s)== w*sm-ImplicitSourceTerm(var=s,coeff=w)-ImplicitSourceTerm(var=s,coeff=un*z1*bn1)-ImplicitSourceTerm(var=s,coeff=ur*z1*br1)-ImplicitSourceTerm(var=s,coeff=un*z2*bn2)-ImplicitSourceTerm(var=s,coeff=ur*z2*br2))
                                       
                                #####################################
                                ######### initial values ############
                                #####################################    
                                
                                br0 = 0.
                                bn0 = 1.
                                bm0 = 1.
                                bm_lim = 5. # unclear how to choose this
                                s0 = 10.
                                i0= 0.
                                r = 6. # radius of initial plug: 3 mm = 6 * delta_x
                                
                                # biomass is placed in a circle with radius r around starting point
                                
                                br1.setValue(br0,where=(x - L/3.)**2 + (y - L/2)**2 < r**2)
                                bn1.setValue(bn0,where=(x - L/3.)**2 + (y - L/2)**2 < r**2)
                                bm1.setValue(bm0,where=(x - L/3.)**2 + (y - L/2)**2 < r**2)
                                
                                br2.setValue(br0,where=(x - (L-L/3.))**2 + (y - L/2)**2 < r**2)
                                bn2.setValue(bn0,where=(x - (L-L/3.))**2 + (y - L/2)**2 < r**2)
                                bm2.setValue(bm0,where=(x - (L-L/3.))**2 + (y - L/2)**2 < r**2)
                                
                                
                                s.setValue(s0)
            
                                # pi = bm/(bn+br), i.e. the ratio of mobile to immobile (= noninsukated and rigid) biomass
                                # using list comprehension to make sure this is done for all values of pi and also to avoid division by zero
                                
                                pi1.setValue([x/y if y else 0. for x,y in zip(bm1,(bn1+br1))],where=(bn1+br1>0)) # /(bn1+bi1)
                                pi2.setValue([x/y if y else 0. for x,y in zip(bm2,(bn2+br2))],where=(bn2+br2>0)) # /(bn1+bi1)
                                
                                Dm1.setValue(Dn1,where=bm1<bm_lim)
                                Dm2.setValue(Dn2,where=bm2<bm_lim)
                                
                                
                                                  
                                
                                #####################################
                                ###### boundary conditions ##########
                                #####################################    
                                
                                valueTopLeft = 0
                                valueBottomRight = 0
                                X, Y = mesh.faceCenters
                                facesTopLeft = ((mesh.facesLeft & (Y > L / 2)) | (mesh.facesTop & (X < L / 2)))
                                facesBottomRight = ((mesh.facesRight & (Y < L / 2)) | (mesh.facesBottom & (X > L / 2)))
                                
                                br1.constrain(valueTopLeft, facesTopLeft)
                                br1.constrain(valueBottomRight, facesBottomRight)
                                bn1.constrain(valueTopLeft, facesTopLeft)
                                bn1.constrain(valueBottomRight, facesBottomRight)
                                bm1.constrain(valueTopLeft, facesTopLeft)
                                bm1.constrain(valueBottomRight, facesBottomRight)
                                
                                br2.constrain(valueTopLeft, facesTopLeft)
                                br2.constrain(valueBottomRight, facesBottomRight)
                                bn2.constrain(valueTopLeft, facesTopLeft)
                                bn2.constrain(valueBottomRight, facesBottomRight)
                                bm2.constrain(valueTopLeft, facesTopLeft)
                                bm2.constrain(valueBottomRight, facesBottomRight)
                                
                                z1.constrain(valueTopLeft, facesTopLeft)
                                z1.constrain(valueBottomRight, facesBottomRight)
                                z2.constrain(valueTopLeft, facesTopLeft)
                                z2.constrain(valueBottomRight, facesBottomRight)
                                
                                s.constrain(valueTopLeft, facesTopLeft)
                                s.constrain(valueBottomRight, facesBottomRight)
                                
                                #-------------------------------------
                                
                                ########################
                                ###### simulation ######
                                ########################
                                
                                
                                eq = eqbr1 & eqbn1 & eqbm1 & eqi1 & eqhsp1 & eqbr2 & eqbn2 & eqbm2 & eqi2 & eqhsp2 & eqs
                                
                                #vi = Viewer(vars=(bn1+br1+bm1+bn2+br2+bm2),FIPY_VIEWER="matplotlib")
                                    
                                ################################
                                ##### start simulation loop ####
                                ################################
                        
                                for loopvar in range(0,101):
                                    
                                    t = loopvar / 10.
                                    
                                    #vi = Viewer(vars=(bn1+br1+bm1+bn2+br2+bm2),FIPY_VIEWER="matplotlib")
                                    
                                    # 1t = 1d -> 1h = 0.04t
                                    
                                    print(t)                       
            
                                    
                                    ##########################
                                    ###### heat effects ######
                                    ##########################
                                    
                                    # turn heat stress on/off:
                                    if t == stressbegin:
                                        temperature.setValue(stresstemp)
                                        
                                        if temperature >= 45.:
                                            z1.setValue(0.,where=(br1+bn1+bm1)>=cutoff)
                                            z2.setValue(0.,where=(br2+bn2+bm2)>=cutoff)
                                            br1.setValue(0.,where=z1==1.)
                                            br2.setValue(0.,where=z2==1.)
                                            bn1.setValue(0.,where=z1==1.)
                                            bn2.setValue(0.,where=z2==1.)
                                            bm1.setValue(0.,where=z1==1.)
                                            bm2.setValue(0.,where=z2==1.)
                                            #z1.setValue(0.)
                                            #z2.setValue(0.)
                                            
                                            h1.setValue((hm1*((temperature-22)/(45-22))),where=br1+bn1+bm1>=cutoff)
                                            h2.setValue((hm2*((temperature-22)/(45-22))),where=br2+bn2+bm2>=cutoff)
                                            
                                    if t == stressbegin+0.1: # ca. 2h later
                                        temperature.setValue(22.)
                                    
                                    h1.setValue(0.,where=((hsp1>hsplim*(br1+bn1+bm1))))
                                    h2.setValue(0.,where=((hsp2>hsplim*(br2+bn2+bm2))))
                                    z1.setValue(1.,where=((hsp1>=hsplim*(br1+bn1+bm1))))
                                    z2.setValue(1.,where=((hsp2>=hsplim*(br2+bn2+bm2))))
                                    
                                
                                    
                                    
                                    #################################
                                    ###### competition effects ######
                                    #################################
                                    
                                    # we do not include autophagy and thus do not include effects of i on alpha/beta
                                    #alpha1.setValue(0.01,where=i2>psi_1)
                                    #alpha2.setValue(0.01,where=i1>psi_2)
                                    #alpha1.setValue(alpha_1,where=i2<=psi_1)
                                    #alpha2.setValue(alpha_2,where=i1<=psi_2)
                                    
                                    #beta1.setValue(0.9,where=i2>psi_1)
                                    #beta2.setValue(0.9,where=i1>psi_2)
                                    #beta1.setValue(beta_1,where=i2<=psi_1)
                                    #beta2.setValue(beta_2,where=i1<=psi_2)
                                    
                                    Dn1.setValue(0.,where=i2>psi_1)
                                    Dn1.setValue(Dn_1,where=i2<=psi_1)
                                    Dn2.setValue(0.,where=i1>psi_2)
                                    Dn2.setValue(Dn_2,where=i1<=psi_2)
                                    
                                
                                    # update pi and Dm
                                    
                                    pi1.setValue([x/y if y else 0. for x,y in zip(bm1,(bn1+br1))],where=(bn1+br1>0)) # /(bn1+bi1)
                                    pi2.setValue([x/y if y else 0. for x,y in zip(bm2,(bn2+br2))],where=(bn2+br2>0)) # /(bn1+bi1)
                                    pi1.setValue(0.,where=(bn1+br1<=0))
                                    pi1.setValue(0.,where=(bm1<=0))
                                    pi2.setValue(0.,where=(bn2+br2<=0))
                                    pi2.setValue(0.,where=(bm2<=0))
                                    
                                    Dm1.setValue(Dn1,where=bm1<bm_lim)
                                    Dm1.setValue(Dn1*10e-7,where=bm1>bm_lim)
                                    Dm2.setValue(Dn2,where=bm2<bm_lim)
                                    Dm2.setValue(Dn2*10e-7,where=bm2>bm_lim)
                                    
                                    
                                        #TSVViewer(vars=(bn1,bi1,n1,pi1,Dn1,z1,hsp1,bn2,bi2,n2,pi2,Dn2,z2,hsp2,s)).plot(filename="dual_hsp_t{0}.txt".format(t))
                                    #vi = Viewer(vars=(n1),FIPY_VIEWER="matplotlib")
                                    
                                    # update all PDEs
                                    br1.updateOld()
                                    br2.updateOld()
                                    bn1.updateOld()
                                    bn2.updateOld()
                                    bm1.updateOld()
                                    bm2.updateOld()
                                    i1.updateOld()
                                    i2.updateOld()
                                    hsp1.updateOld()
                                    hsp2.updateOld()
                                    s.updateOld()
                                    
                                    if (t == 10): #or ( t > 2 and t < 3.5):
                                        # xx is a temporary dataframe that lists all variables and then appends them to DF
                                        xx = CustomViewer(vars=(bn1,br1,bm1,pi1,Dn1,z1,hsp1,i1,h1,bn2,br2,bm2,pi2,Dn2,z2,hsp2,i2,h2,s)).makedf(time=t)
                                        xx.insert(2,"Db1",i)
                                        xx.insert(2,"Db2",j)
                                        xx.insert(2,"hm1",k)
                                        xx.insert(2,"hm2",l)
                                        xx.insert(2,"Di1",m)
                                        xx.insert(2,"Di2",q)
                                        xx.insert(2,"stresstemp",te)
                                        DF = DF.append(xx)
                                    
                                        # show plot of the biomass dynamics
                                    #vi = Viewer(vars=(bm1+br1+bn1+bm2+br2+bn2),FIPY_VIEWER="matplotlib")
                                    
                                    
                                    eq.solve(dt=0.1)
                                    
                                    
                                #y.to_csv(("noheat_omega%d.txt" % (i*10)),index=False)
                                
                            
                    
            # save final dataframe 
        
DF.to_csv("alloutcomes.txt",index=False)
    


