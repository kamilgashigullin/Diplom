import numpy as np
from constants import *
def fsolve(FF,ne,ni,omega_relax = 1.8):
    FF1 = FF.copy()
   # FF_glob = np.reshape(FF,N*M)
    f =   np.zeros((N,M)) #e/eps0*(ni-ne)
   # f_liq = 0
    d = np.ones((N*M)) #формируем диагональ матрицы СЛАУ
    l1 = np.zeros((N*M))
    l1_ = np.zeros((N*M))
    l2 = np.zeros((N*M))
    l2_ = np.zeros((N*M))

    '''
    l1 = np.zeros((N*M-1)) #верхняя блихняя диагональ
    l1_ = np.zeros((N*M-1)) #нижняя ближняя диагональ
    l2 = np.zeros((N*M-N)) #верхняя дальняя диагональ
    l2_ = np.zeros((N*M-N)) #нижхняя дальняя диагональ
    '''
    for j in range(1,M-1):
        for i in range(N-1):
            k = i+j*N
            d[k] = 2/hr**2 + 2/hz**2 # главная диагональ (r[i]+hr/2)/(r[i]*hr**2)
            l1[k] = -(r[i]+hr/2)/(r[i]*hr**2)
            l1_[k-1] = -(r[i]-hr/2)/(r[i]*hr**2)
            l2[k] = -1/hz**2
            l2_[k-N] = -1/hz**2
            print("Summa norma =     ", l1[k]+l1_[k-1]+l2[k]+l2_[k-N]+d[k], 'k =  ', k )


    tol = (1e-2)
    rn_finish = (1e32)

    kiter = 0
    N_iter = 5

    for p in range(N_iter):
        norm_dif = 0
        rn = 0
        print('rn = ',float(rn))

        print("Расчет")
      
        resid = np.zeros((N*M))
        corrector = np.zeros((N*M))

        
        for j in range(1,M-1):
            # Граничные условия второго рода слева  
            k = j*N
           #print('Слагаемые  для вычисления невязки на левой границе')
            #print("d[k] = ",d[k],'FF[0,j] =  ',FF[0,j],'l1[k] = ',l1[k],'FF[1,j] = ',FF[1,j], 'l2[k] = ',l2[k],'FF[0,j+1] = ',FF[0,j+1],'f[0,j] = ',f[0,j], 'k =  ',k )
            rr = d[k]*FF[0,j] + l1[k]*FF[1,j] + l2[k]*FF[0,j+1]+ l2_[k-N]*FF[0,j-1] - f[0,j]          
            # rr = (FF[0,j]*(2/hr**2) - FF[1,j]*(2/hr**2))\
            # - 1/hz**2 * FF[0,j-1] + 2/hz**2 * FF[0,j] - 1/hz**2 * FF[0,j+1]\
            # - f[0,j]
            resid[j*N] = rr
            rn = rn+rr**2
            print('rr невязка левые граница', rr )
            FF1[0,j] = FF[0,j] - (omega_relax*rr)/d[j*N]
            norm_dif = max(norm_dif,abs((omega_relax*rr)/(d[j*N])))
            
            # Обходим внутренние узлы
            for i in range(1,N-1):
                k = i+j*N
                print("Слагаемые для расчета невязки по внутренним узлам")
                print("d[k] = ",d[k],'FF[i,j] =  ',FF[i,j],'l1[k] = ',l1[k],'l1_[k-1] = ',l1_[k-1],'FF[i-1,j] = ',FF[i-1,j], 'FF[i+1,j] = ',FF[i+1,j], 'l2[k] = ',l2[k],'FF[i,j+1] = ',FF[i,j+1],'l2_[k-N] = ',l2_[k-N],'FF[i,j-1] = ',FF[i,j-1],'f[i,j] = ',f[i,j], 'k =  ',k )
                rr = l2_[k-N]*FF[i,j-1] + l1_[k-1]*FF[i-1,j] + d[k]*FF[i,j] + l1[k]*FF[i+1,j] + l2[k]*FF[i,j+1] - f[i,j]
                # rr = -(FF[i-1,j]*(r[i]- hr/2) - 2*FF[i,j]*r[i] + FF[i+1,j]*(r[i]+hr/2)) / (hr**2 * r[i]) -(FF[i,j-1] - 2*FF[i,j] + FF[i,j+1])/hz**2 - f[i,j]
                #rr_glob = d[i]/FF_glob[i] + l1_[i]*FF_glob[i-N-1] + 
                corrector[i+j*N] = rr - omega_relax*(l1[i+j*N-1]*corrector[i+j*N -1] 
                                                     - l2[i+(j-1)*N]*corrector[i+(j-1)*N])/d[i+j*N]
              
                resid[i+j*N] = rr
                rn = rn+rr**2
                print ('i = ',i,'j = ', j, 'rr = ',float(rr),'rn = ',float(rn))
                FF1[i,j] = FF[i,j] - omega_relax*corrector[i+j*N]
                norm_dif = max(norm_dif,abs(omega_relax*corrector[i+j*N]))

            # справа
          #  k = (j+1)*N-1
           # rr  = l2_[k-N]*FF[N-1,j-1] + l1_[k-1]*FF[N-2,j] + d[k]*FF[N-1,j]  - f[N-1,j]
            # rr = (FF[N-1,j] - FF[N-2,j])*2 / hr**2 - f[N-1,j]\
                #- 1/hz**2 * (FF[N-1,j-1] - 2*FF[N-1,j] + FF[N-1,j+1])
          #  corrector[(j+1)*N-1] = rr - omega_relax*(l1[(j+1)*N-2]*corrector[(j+1)*N-2] 
                                              #      - l2[j*N-1]*corrector[j*N-1])/d[(j+1)*N-1]
            #print('невязка правой границы',float(rr),'\n')
           # resid[(j+1)*N-1] = rr
            #rn = rn+rr**2
            #FF1[N-1,j] = FF[N-1,j] - (omega_relax*corrector[(j+1)*N-1])
           # norm_dif = max(norm_dif,abs(omega_relax*corrector[(j+1)*N-1]))

        FF=FF1.copy()
        rn = rn*hr*hz
        kiter = p
        rn_finish = rn
        #print("итерация номер ", k, " невязка = ", rn," погрешность = ", err_c)
        if rn_finish < tol**2:
            break
      
    print("число итераций =  ", kiter+1)
    print("невязка приближенного решения",rn)
    print("относительная разность, norm_dif = ", norm_dif)
    return FF,rn
