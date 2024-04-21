import matplotlib.animation
import numpy as np
import matplotlib.pyplot as plt


class fourBar():

    def __init__(self):
        pass

    def pos_solver(self, a, b, c, d,
                            ta, tb, mode=0):

        x = a*np.cos(ta) + b*np.cos(tb)
        y = a*np.sin(ta) + b*np.sin(tb)

        R = np.sqrt(x**2 + y**2)

        Ac = (c**2 - d**2 + R**2)/(2*c)
        Ad = (d**2 - c**2 + R**2)/(2*d)

        tc = 2*np.arctan2(- y + ((-1)**mode)*np.sqrt(R**2 - Ac**2), -x-Ac)
        td = 2*np.arctan2(- y - ((-1)**mode)*np.sqrt(R**2 - Ad**2), -x-Ad)

        return tc, td

    def vel_solver(self, a, b, c, d,
                            ta, tb, tc, td,
                            wa, wb):

        A = np.array([[c*np.cos(tc), d*np.cos(td)],
                      [c*np.sin(tc), d*np.sin(td)]])

        y = np.array([a*wa*np.cos(ta)+b*wb*np.cos(tb),
                      a*wa*np.sin(ta)+b*wb*np.sin(tb)])

        wc, wd = np.linalg.solve(A, y.T)
        return wc, wd

    def acc_solver(self, a, b, c, d,
                            ta, tb, tc, td,
                            wa, wb, wc, wd,
                            aa, ab):

        A = np.array([[c*np.cos(tc), d*np.cos(td)],
                      [c*np.sin(tc), d*np.sin(td)]])

        y = np.array([(a*aa*np.cos(ta) - a*(wa**2)*np.sin(ta) + c*(wc**2)*np.sin(tc)
                       + b*ab*np.cos(tb) - b*(wb**2)*np.sin(tb) + d*(wd**2)*np.sin(td)),
                      (a*aa*np.sin(ta) + a*(wa**2)*np.cos(ta) - c*(wc**2)*np.cos(tc)
                     + b*ab*np.sin(tb) + b*(wb**2)*np.cos(tb) - d*(wd**2)*np.cos(td))])

        ac, ad = np.linalg.solve(A, y.T)
        return ac, ad

    def solver(self, a, b, c, d,
                        ta, tb,
                        wa, wb,
                        aa, ab,
                        mode=0):

        tc, td = self.pos_solver(a, b, c, d,
                                 ta, tb, mode=mode)

        wc, wd = self.vel_solver(a, b, c, d,
                                 ta, tb, tc, td,
                                 wa, wb)

        ac, ad = self.acc_solver(a, b, c, d,
                                 ta, tb, tc, td,
                                 wa, wb, wc, wd,
                                 aa, ab)

        return tc, td, wc, wd, ac, ad

    def plot(self, a, b, c, d,
                      ta, tb, tc, td,
                      ox=0, oy=0):

        plt.grid("on")
        plt.plot([ox+0, ox+a*np.cos(ta)],
                 [oy+0, oy+a*np.sin(ta)])

        plt.plot([ox+a*np.cos(ta), ox+a*np.cos(ta)+b*np.cos(tb)],
                 [oy+a*np.sin(ta), oy+a*np.sin(ta)+b*np.sin(tb)])

        plt.plot([ox+0, ox+d*np.cos(td)],
                 [oy+0, oy+d*np.sin(td)])

        plt.plot([ox+d*np.cos(td), ox+d*np.cos(td)+c*np.cos(tc)],
                 [oy+d*np.sin(td), oy+d*np.sin(td)+c*np.sin(tc)])


class Point():
    def __init__(self, a=0, ta=0, wa=0, aa=0):

        self.pos = np.array([a*np.cos(ta), a*np.sin(ta)])
        self.vel = np.array([-a*wa*np.sin(ta), a*wa*np.cos(ta)])
        self.acc = np.array([-a*(wa**2)*np.cos(ta) - a*aa*np.sin(ta),
                             -a*(wa**2)*np.sin(ta) + a*(aa)*np.cos(ta)])

    def __add__(self, b):
        c = Point()
        c.pos = self.pos + b.pos
        c.vel = self.vel + b.vel
        c.acc = self.acc + b.acc
        return c

    def __sub__(self, b):
        c = Point()
        c.pos = self.pos - b.pos
        c.vel = self.vel - b.vel
        c.acc = self.acc - b.acc
        return c

    def __neg__(self):
        c = Point()
        c.pos = - self.pos
        c.vel = - self.vel
        c.acc = - self.acc
        return c
    
    def __mul__(self, n):
        c = Point()
        c.pos = n*self.pos
        c.vel = n*self.vel
        c.acc = n*self.acc
        return c

    def __truediv__(self, n):
        c = Point()
        c.pos = self.pos/n
        c.vel = self.vel/n
        c.acc = self.acc/n
        return c       

    def __repr__(self):
        str = f"__name__\n"
        str += f"Position    : {self.pos}\n"
        str += f"Velocity    : {self.vel}\n"
        str += f"Acceleration: {self.acc}\n"
        return str


def freudenstein_method(phi, psi, r1):

    # psi = tht4, phi = tht2
    
    # k0 - k2 c_tht4 + k4 c_tht2 = c_(tht4-tht2)
    
    A = np.array([[1, -np.cos(psi[0]), np.cos(phi[0])],
                  [1, -np.cos(psi[1]), np.cos(phi[1])],
                  [1, -np.cos(psi[2]), np.cos(phi[2])]])
    
    b = np.array([np.cos(phi[0]-psi[0]),
                  np.cos(phi[1]-psi[1]),
                  np.cos(phi[2]-psi[2])])
    
    k0, k2, k4 = np.linalg.solve(A, b)
    
    
    r2 = r1/k2
    r4 = r1/k4
    r3 = np.sqrt(r1**2 + r2**2 + r4**2 - 2*r2*r4*k0)

    return r1, r2, r3, r4

def d_fn(g1, g2, d1, d2):

    # | g1 d1 |
    # | g2 d2 |
    
    D = np.array([   (np.cos(g1)-1)*d2[0] - d2[1]*np.sin(g1)
                   -((np.cos(g2)-1)*d1[0] - d1[1]*np.sin(g2)),
                     (np.cos(g1)-1)*d2[1] + d2[0]*np.sin(g1)
                   -((np.cos(g2)-1)*d1[1] + d1[0]*np.sin(g2))])
    return D
    
def burmester_method(gam, delta, phi2, mode=0, debug=False):

    # delta[0] = delta_2 ... delta[2] = delta_4
    # gam[0] = gamma_2 ... gam[2] = gamma_4
    # phi[0] = phi_2 ... phi[2] = phi_4
    
    D2 =  d_fn(gam[1], gam[2], delta[1], delta[2])
    D3 = -d_fn(gam[0], gam[2], delta[0], delta[2])
    D4 =  d_fn(gam[0], gam[1], delta[0], delta[1])
    
    D1 = -D2 -D3 -D4
    D = D1 + D2@np.array([[ np.cos(phi2), np.sin(phi2)],
                          [-np.sin(phi2), np.cos(phi2)]])


    [nD, nD1, nD2, nD3, nD4] = [np.linalg.norm(di) for di in [D, D1, D2, D3, D4]]

    if(((nD+nD3)>nD4) and ((nD4+nD)>nD3) and ((nD3+nD4)>nD)):
        pass
    else:
        if debug:
            print(f"ERROR: no solution for phi3, phi4, given phi2 {phi2}")
        return np.array([np.nan, np.nan]), np.array([np.nan, np.nan]), np.array([np.nan, np.nan, np.nan])
        
    c_tht3 = -(nD3**2 + nD**2 - nD4**2)/(2*nD3*nD)    
    c_tht4 = -(nD4**2 + nD**2 - nD3**2)/(2*nD4*nD)
    
    tht3 = ((-1)**(mode))*np.arccos(c_tht3)
    tht4 = ((-1)**(mode))*np.arccos(c_tht4)

    phi3 = np.arctan2(D[1], D[0]) - np.arctan2(D3[1], D3[0]) + tht3
    phi4 = np.arctan2(D[1], D[0]) - np.arctan2(D4[1], D4[0]) - tht4

    phi = np.array([phi2, phi3, phi4])

    A = np.array([[np.cos(phi2)-1, -np.sin(phi2),   np.cos(gam[0])-1, -np.sin(gam[0])],
                  [np.sin(phi2)  ,  np.cos(phi2)-1, np.sin(gam[0])  ,  np.cos(gam[0])-1],
                  [np.cos(phi3)-1, -np.sin(phi3),   np.cos(gam[1])-1, -np.sin(gam[1])],
                  [np.sin(phi3)  ,  np.cos(phi3)-1, np.sin(gam[1])  ,  np.cos(gam[1])-1]])

    b = np.array([delta[0,0], 
                  delta[0,1],
                  delta[1,0], 
                  delta[1,1]])

    try:
        z1x, z1y, z2x, z2y = np.linalg.solve(A, b)
    except:
        if debug:
            print(f"ERROR: Singular matrix for phi2 {phi2}")
        return np.array([np.nan, np.nan]), np.array([np.nan, np.nan]), np.array([np.nan, np.nan, np.nan])
        

    z1 = np.array([z1x, z1y])
    z2 = np.array([z2x, z2y])

    return z1, z2, phi 

if __name__ == "__main__": 
    pass