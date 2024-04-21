import matplotlib.animation
import numpy as np
import matplotlib.pyplot as plt

from tom_utils import fourBar
from tom_utils import Point


class jansenLinkage():

    def __init__(self, links=None, tn=None):

        self.link = {"b": 41.5, "c": 39.3,
                     "d": 40.1, "e": 55.8, "f": 39.4,
                     "g": 36.7, "h": 65.7, "i": 49.0,
                     "j": 50.0, "k": 61.9, "n": (38.0**2 + 7.8**2)**0.5,
                     "m": 15.0}

        self.tn = np.arctan2(7.8, 38.0)

        if links is not None:
            self.link = links

        if tn is not None:
            self.tn = tn

    def kinematic_analysis(self, tm, wm=0, wn=0, am=0, an=0):

        tn = self.tn

        # eqn 1
        tj, tb, wj, wb, aj, ab = fourBar().solver(self.link["n"], self.link["m"],
                                                  self.link["j"], self.link["b"],
                                                  self.tn, tm,
                                                  wn, wm,
                                                  an, am,
                                                  mode=0)
        # eqn 2
        tk, tc, wk, wc, ak, ac = fourBar().solver(self.link["n"], self.link["m"],
                                                  self.link["k"], self.link["c"],
                                                  self.tn, tm,
                                                  wn, wm,
                                                  an, am,
                                                  mode=1)
        # eqn 4
        td = tb + np.arccos((self.link["b"]**2+self.link["d"]**2 -
                            self.link["e"]**2)/(2*self.link["b"]*self.link["d"])) + np.pi
        wd = wb
        ad = ab

        # eqn 3
        tg, tf, wg, wf, ag, af = fourBar().solver(self.link["d"], self.link["c"],
                                                  self.link["g"], self.link["f"],
                                                  td, tc,
                                                  wd, wc,
                                                  ad, ac,
                                                  mode=1)
        # eqn 5
        ti = tg + np.arccos((self.link["g"]**2+self.link["i"]**2 -
                            self.link["h"]**2)/(2*self.link["g"]*self.link["i"])) - np.pi
        wi = wg
        wh = wg
        ai = ag
        ah = ag

        # eqn 6
        px = self.link["c"]*np.cos(tc) + self.link["i"]*np.cos(ti)
        py = self.link["c"]*np.sin(tc) + self.link["i"]*np.sin(ti)

        pos = {"tb": tb, "tc": tc, "td": td,
               "tf": tf, "tg": tg, "ti": ti,
               "tj": tj, "tk": tk, "tn": self.tn,
               "tm": tm, "px": px, "py": py}

        ang_vel = {"wb": wb, "wc": wc, "wd": wd,
                   "wf": wf, "wg": wg, "wi": wi,
                   "wj": wj, "wk": wk, "wn": wn,
                   "wm": wm}

        ang_acc = {"ab": ab, "ac": ac, "ad": ad,
                   "af": af, "ag": ag, "ai": ai,
                   "aj": aj, "ak": ak, "an": an,
                   "am": am}

        # declare points for kinematic analysis
        vecs = {}
        for ci in ["b", "c", "d", "f", "g", "i", "j", "k", "n", "m"]:
            vecs[ci] = Point(
                self.link[ci], pos[f"t{ci}"], ang_vel[f"w{ci}"], ang_acc[f"a{ci}"])

        points = {"p0": Point(0, 0, 0, 0),      "p3": vecs["n"],
                  "p1": vecs["n"]+vecs["m"], "p2": vecs["b"],
                  "p7": vecs["c"]+vecs["i"], "p5": vecs["c"],
                  "p6": vecs["c"]-vecs["g"], "p4": -vecs["d"]}

        return pos, ang_vel, ang_acc, points

    def plot_jansen_linkage(self, idx):
        plt.cla()

        plt.xlim([-160, 120])
        plt.ylim([-200, 120])
        plt.gca().set_aspect('equal')

        res, _, _, _ = self.kinematic_analysis(2*idx*np.pi/180)

        fourBar().plot(self.link["n"], self.link["m"],
                       self.link["j"], self.link["b"],
                       res["tn"], res["tm"], res["tj"], res["tb"])

        fourBar().plot(self.link["n"], self.link["m"],
                       self.link["k"], self.link["c"],
                       res["tn"], res["tm"], res["tk"], res["tc"])

        fourBar().plot(self.link["d"], self.link["c"],
                       self.link["g"], self.link["f"],
                       res["td"], res["tc"], res["tg"], res["tf"],
                       ox=-(self.link["d"]*np.cos(res["td"])),
                       oy=-(self.link["d"]*np.sin(res["td"])))

        plt.plot([self.link["c"]*np.cos(res["tc"]), res["px"]],
                 [self.link["c"]*np.sin(res["tc"]), res["py"]])

        plt.plot([self.link["c"]*np.cos(res["tc"])-self.link["g"]*np.cos(res["tg"]), res["px"]],
                 [self.link["c"]*np.sin(res["tc"])-self.link["g"]*np.sin(res["tg"]), res["py"]])

        plt.plot([self.link["b"]*np.cos(res["tb"]), -self.link["d"]*np.cos(res["td"])],
                 [self.link["b"]*np.sin(res["tb"]), -self.link["d"]*np.sin(res["td"])])

        plt.plot(res["px"], res["py"], "ro")
        plt.plot(self.couper_curve[:, 0],self.couper_curve[:, 1])

    def display(self, file="../outputs/jansen_linkage.gif", save=True):

        # generate coupler curve
        self.couper_curve  = []
        for i in range(360):
            pos, _, _ , _ = self.kinematic_analysis(tm=i*np.pi/180)
            self.couper_curve.append([pos["px"], pos["py"]])

        self.couper_curve = np.array(self.couper_curve)       

        self.fig, self.ax = plt.subplots(figsize=(10, 5))
        self.animation = matplotlib.animation.FuncAnimation(self.fig,
                                                            self.plot_jansen_linkage,
                                                            frames=180)

        writer = matplotlib.animation.PillowWriter(fps=72)
        if save:
            self.animation.save(file, writer=writer)
        
        return self.animation


if __name__ == "__main__":
    pass
