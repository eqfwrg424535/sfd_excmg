# script for generating the input files, using DTM1.O as an example
# author: Jinxuan Wang (06Dec2024)
import numpy as np

#function for generating inputf:
def gen_inputfor(modelname):
     with open(modelname+'.fwd','w') as input:  
         print('The DTM1.0 model generated with Python3',file = input)
         level = 4 # number of grid layers, including the initial two coarse grids
         eps = 1e-10 # tolerance the finest grid level
         solver = 1 # solver type: 1--BiCGStab; 2--TFQMR; 3--GPBiCG
         out_type = 1 # 1--impeadance(real and imag); 2--apparent resistivity and phase 
        # precond = 1 # preconditioner type: 1--SSOR; 2--ILU(0); 3--Jacobi
         print('%d %e %d %d' %(level, eps, solver, out_type),file =input)
        # the number of cells in each direction, the air layers, and the indicator for topography
         nx = 128; ny = 128; nz = 128; nair = 64; topo_indic = 0  
         print('%d %d %d %d %d' %(nx,ny,nz,nair,topo_indic), file = input) 
         a0 = np.zeros(16);b0=np.zeros(16);c0=np.zeros(16)
         a0[4:12] = 5e3; a0[:4] = [2.92e4,1.62e4,9e3,5e3]
         a0[12:16] = a0[:4][::-1]
         b0 = a0; b0[3] = 5e3; b0[12]=5e3
         c0[:10] = 5e3
         for i in range(6):
             c0[10+i] = 5e3*np.power(1.8,i+1)
         ax = np.zeros(nx)
         by = np.zeros(ny)
         cz = np.zeros(nz+nair)
         ext_h = 1.15
         ext_v = 1.08
         for i in range(16):
             ax[i*8:(i+1)*8] = a0[i]/8
             by[i*8:(i+1)*8] = b0[i]/8
             cz[nair+i*8:nair+(i+1)*8] = c0[i]/8
         cz[:nair] = cz[128:][::-1]
         for i in range(nx):
             if i !=nx-1:
                 print('%.4e' %ax[i], end=' ', file = input)
             else:
                 print('%.4e' %ax[i], end='\n', file = input)
         for i in range(ny):
             if i !=ny-1:
                 print('%.4e' %by[i], end=' ', file = input)
             else:
                 print('%.4e' %by[i], end='\n', file = input)
         for i in range(nair,nz+nair):
             if i !=nz+nair-1:
                 print('%.4e' %cz[i], end=' ', file = input)
             else:
                 print('%.4e' %cz[i], end='\n', file = input)
 
         # the extension factor of both horizontal and vertical directions (no bigger than 2)
         print('%f %f' %(ext_h,ext_v),file = input)
         # the array indicating the range of uniform cells(considering air layers)
         #abu = np.array([65,192,65,192,57,192])
         abu = np.array([25,104,25,104,65,144])
         for i in range(6):
             if i !=5:
                 print('%d' %(abu[i]), end=' ', file = input)
             else:
                 print('%d' %(abu[i]), end='\n', file = input)
         #print('%d %d %d %d %d %d' %(abu),file = input)
         ncell = nx*ny*nz
         conduc = np.zeros([ncell,6]) # column 1~3 are for main-axis conductivities, while 4~6 are for off-digonal
         conduc[:,:2] = 1e-2 # conductivity of half-space
         # list that defines the locations of the anomalous bodies, only for simple synthetic models
        # abr = [88,168,123,133,10,50, 98,128,83,133,50,60, 128,158,123,173,50,100]
         abr = [32,96,60,68,8,32, 40,64,60,100,32,40, 64,88,28,68,32,80]
         anb_conduc = [0.1,1,1e-4]
         # indicator of the conductivity type (1--isotropic; 2--axial anisotropic; 3--arbitary anisotropic)
         # to save the time of file reading
         conduc_type = 1
         print('%d' %(conduc_type), file = input)
         # deal with the conductivities 
         for i in range(int(len(abr)/6)):
              xstart = abr[6*i]; xend = abr[6*i+1]
              ystart = abr[6*i+2]; yend = abr[6*i+3]
              zstart = abr[6*i+4]; zend = abr[6*i+5]
              for l in range(zstart,zend):
                   for m in range(ystart,yend):
                        for n in range(xstart,xend):
                            conduc[l*nx*ny+m*nx+n,:3] = anb_conduc[i]
            #  for j in range(6):
         for i in range(np.size(conduc,0)):
            if conduc_type == 1:
                  print('%.4e' %(conduc[i,0]), file = input)
            elif conduc_type == 2:
                  print('%.4e %.4e %.4e' %(conduc[i,0],conduc[i,1],conduc[i,2]), file = input)
            else:
                for j in range(6):
                    if j !=5:
                        print('%.4e' %(conduc[i,j]), end=' ', file = input)
                    else:
                        print('%.4e' %(conduc[i,j]), end='\n', file = input)
              #print('%e %e %e %e %e %e' %(conduc[i,:]),file = input)


def gen_datfile(modelname):
     nf = 4
     freq = [0.001, 0.01, 0.1, 1]
     nsite = 17 #4 lines
     xsite = [-37.5, -32.5, -27.5, -22.5, -17.5, -12.5, -7.5, -2.5, 0.0, \
                2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5]
     #xsite = xsite*1000
     ysite = 0; zsite = 0
     ysite2 = 15000; ysite3 = -15000
     with open(modelname+'.dat','w') as input:  
         print('%d' %nf, file = input)
         for i in range(nf):
              print('%f' %freq[i], file = input)
         print('%d' %nsite, file = input)
         for i in range(nsite): 
              print('%e %e %e' %(xsite[i]*1000,ysite,zsite), file = input)
         for i in range(nsite): 
              print('%e %e %e' %(ysite3,xsite[i]*1000,zsite), file = input)
         for i in range(nsite): 
              print('%e %e %e' %(ysite,xsite[i]*1000,zsite), file = input)
         for i in range(nsite): 
              print('%e %e %e' %(ysite2,xsite[i]*1000,zsite), file = input)

# function for calculating the conductivity tensor using Euler matrixes
# input: conduc is an 1*6 array with (sigma_x,sigma_y,sigma_z,alpha_s,alpha_d,alpha_l)
# output: sigma is the 3*3 sigma tensor
def conduc_cal(conduc):
    main_con = np.diag(np.array(conduc[:3]))
    alphas = conduc[3] # strike angle
    alphad = conduc[4] # dip angle
    alphal = conduc[5] # slant angle 
    matrixs = np.array([[np.cos(alphas), -np.sin(alphas), 0],
                   [np.sin(alphas), np.cos(alphas), 0],
                   [0, 0, 1]])
    matrixd = np.array([[1, 0, 0],
                   [0, np.cos(alphad), -np.sin(alphad)],
                   [0, np.sin(alphad), np.cos(alphad)]])
    matrixl = np.array([[np.cos(alphal), -np.sin(alphal), 0],
                   [np.sin(alphal), np.cos(alphal), 0],
                   [0, 0, 1]])
    dtmp = np.matmul(matrixs,np.matmul(matrixd,matrixl))
    sigma = np.matmul(dtmp,np.matmul(main_con,np.transpose(dtmp)))

    return sigma


if __name__ == "__main__":
     modelname = "dtm1"
     gen_inputfor(modelname)
     gen_datfile(modelname)