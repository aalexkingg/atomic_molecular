from vpython import *
#import base.atomic as atomic

# constants
WIDTH = 800
HEIGHT = 400

k = 1.4E-23  # Boltzmann constant


class GasAnimation:
    def __init__(self, atom=None):
        """
        :param atom: Atom Element
        """
        self.deltav = 100  # binning for v histogram
        self.nhisto = 0  # number of histogram snapshots to average

        self.animation = canvas(width=WIDTH / 3, height=HEIGHT, align='left')
        self.animation.userpan = True

        # self.particle_count_slider = slider(bind=set_count, min=1, max=300, align='right')

        self.num_atoms = 200  # change this to have more or fewer atoms

        # Typical values
        self.L = 1  # container is a cube L on a side
        self.gray = color.gray(0.7)  # color of edges of container
        self.mass = 4E-3 / 6E23  # helium mass
        self.Ratom = 0.03  # wildly exaggerated size of helium atom

        self.T = 300  # around room temperature
        self.dt = 1E-5

        self.atoms = []
        self.p = []
        self.apos = []
        self.p_avg = sqrt(2 * self.mass * 1.5 * k * self.T)  # average kinetic energy p**2/(2mass) = (3/2)kT

        self.nhisto = int(4500 / self.deltav)
        self.histo = []

        self.draw_box()
        self.draw_speed_graph()

    def start(self):
        """

        :return:
        """
        self.create_particles()

        for i in range(self.nhisto):
            self.histo.append(0.0)

        self.histo[self.barx(self.p_avg / self.mass)] = self.num_atoms

        accum = []
        for i in range(int(3000 / self.deltav)):
            accum.append([self.deltav * (i + .5), 0])

        vdist = gvbars(color=color.red, delta=self.deltav)

        while True:
            # Check for input event from sliders/inputs



            rate(300)

            # Accumulate and average histogram snapshots
            for i in range(len(accum)):
                accum[i][1] = (self.nhisto * accum[i][1] + self.histo[i]) / (self.nhisto + 1)

            if self.nhisto % 10 == 0:
                vdist.data = accum
            self.nhisto += 1

            # Update all positions
            for i in range(self.num_atoms):
                self.atoms[i].pos = self.apos[i] = self.apos[i] + (self.p[i] / self.mass) * self.dt

            # Check for collisions
            # If any collisions took place, update momenta of the two atoms
            for ij in self.checkCollisions():
                i = ij[0]
                j = ij[1]
                ptot = self.p[i] + self.p[j]
                posi = self.apos[i]
                posj = self.apos[j]
                vi = self.p[i] / self.mass
                vj = self.p[j] / self.mass
                vrel = vj - vi
                a = vrel.mag2
                if a == 0:
                    continue  # exactly same velocities

                rrel = posi - posj
                if rrel.mag > self.Ratom:
                    continue  # one atom went all the way through another

                # theta is the angle between vrel and rrel:
                dx = dot(rrel, vrel.hat)  # rrel.mag*cos(theta)
                dy = cross(rrel, vrel.hat).mag  # rrel.mag*sin(theta)
                # alpha is the angle of the triangle composed of rrel, path of atom j, and a line
                #   from the center of atom i to the center of atom j where atome j hits atom i:
                alpha = asin(dy / (2 * self.Ratom))
                d = (2 * self.Ratom) * cos(alpha) - dx  # distance traveled into the atom from first contact
                deltat = d / vrel.mag  # time spent moving from first contact to position inside atom

                posi = posi - vi * deltat  # back up to contact configuration
                posj = posj - vj * deltat
                mtot = 2 * self.mass
                pcmi = self.p[i] - ptot * self.mass / mtot  # transform momenta to cm frame
                pcmj = self.p[j] - ptot * self.mass / mtot
                rrel = norm(rrel)
                pcmi = pcmi - 2 * pcmi.dot(rrel) * rrel  # bounce in cm frame
                pcmj = pcmj - 2 * pcmj.dot(rrel) * rrel
                self.p[i] = pcmi + ptot * self.mass / mtot  # transform momenta back to lab frame
                self.p[j] = pcmj + ptot * self.mass / mtot
                self.apos[i] = posi + (self.p[i] / self.mass) * deltat  # move forward deltat in time
                self.apos[j] = posj + (self.p[j] / self.mass) * deltat
                self.interchange(vi.mag, self.p[i].mag / self.mass)
                self.interchange(vj.mag, self.p[j].mag / self.mass)

            for i in range(self.num_atoms):
                loc = self.apos[i]
                if abs(loc.x) > self.L / 2:
                    if loc.x < 0:
                        self.p[i].x = abs(self.p[i].x)
                    else:
                        self.p[i].x = -abs(self.p[i].x)

                if abs(loc.y) > self.L / 2:
                    if loc.y < 0:
                        self.p[i].y = abs(self.p[i].y)
                    else:
                        self.p[i].y = -abs(self.p[i].y)

                if abs(loc.z) > self.L / 2:
                    if loc.z < 0:
                        self.p[i].z = abs(self.p[i].z)
                    else:
                        self.p[i].z = -abs(self.p[i].z)

    def restart(self):
        """

        :return:
        """
        ...
        self.start()

    def draw_box(self):
        """

        :return:
        """
        d = self.L / 2 + self.Ratom
        r = 0.005
        boxbottom = curve(color=self.gray, radius=r)
        boxbottom.append([vector(-d, -d, -d), vector(-d, -d, d), vector(d, -d, d), vector(d, -d, -d), vector(-d, -d, -d)])
        boxtop = curve(color=self.gray, radius=r)
        boxtop.append([vector(-d, d, -d), vector(-d, d, d), vector(d, d, d), vector(d, d, -d), vector(-d, d, -d)])
        vert1 = curve(color=self.gray, radius=r)
        vert2 = curve(color=self.gray, radius=r)
        vert3 = curve(color=self.gray, radius=r)
        vert4 = curve(color=self.gray, radius=r)
        vert1.append([vector(-d, -d, -d), vector(-d, d, -d)])
        vert2.append([vector(-d, -d, d), vector(-d, d, d)])
        vert3.append([vector(d, -d, d), vector(d, d, d)])
        vert4.append([vector(d, -d, -d), vector(d, d, -d)])

    def create_particles(self):
        """

        :return:
        """
        for i in range(self.num_atoms):
            x = self.L * random() - self.L / 2
            y = self.L * random() - self.L / 2
            z = self.L * random() - self.L / 2
            if i == 0:
                self.atoms.append(
                    sphere(
                        pos=vector(x, y, z),
                        radius=self.Ratom / 2,
                        color=color.cyan,
                        make_trail=False,
                        retain=100,
                        trail_radius=0.3 * self.Ratom
                    )
                )
            else:
                self.atoms.append(sphere(pos=vector(x, y, z), radius=self.Ratom / 2, color=self.gray))
            self.apos.append(vec(x, y, z))
            theta = pi * random()
            phi = 2 * pi * random()
            px = self.p_avg * sin(theta) * cos(phi)
            py = self.p_avg * sin(theta) * sin(phi)
            pz = self.p_avg * cos(theta)
            self.p.append(vector(px, py, pz))

    def draw_speed_graph(self):
        """

        :return:
        """
        speed_graph = graph(
            width=WIDTH / 3,
            height=HEIGHT / 2,
            xmax=3000,
            align='right',
            xtitle='speed, m/s',
            ytitle='Number of atoms',
            ymax=100
        )
        theory = gcurve(color=color.blue, width=2)

        dv = 10
        for v in range(0, 3001 + dv, dv):  # theoretical prediction
            theory.plot(v, (self.deltav / dv) * self.num_atoms * 4 * pi * ((self.mass / (2 * pi * k * self.T)) ** 1.5) * exp(
                -0.5 * self.mass * (v ** 2) / (k * self.T)) * (v ** 2) * dv)

    def checkCollisions(self):
        """

        :return:
        """
        hitlist = []
        r2 = 2 * self.Ratom
        r2 *= r2
        for i in range(self.num_atoms):
            ai = self.apos[i]
            for j in range(i):
                aj = self.apos[j]
                dr = ai - aj
                if mag2(dr) < r2: hitlist.append([i, j])
        return hitlist




    def barx(self, v):
        """

        :param v:
        :return:
        """
        return int(v / self.deltav)  # index into bars array


    def set_count(self, evt):
        """

        :param evt:
        :return:
        """
        self.num_atoms = floor(evt.value)


    def interchange(self, v1, v2):  # remove from v1 bar, add to v2 bar
        """

        :param v1:
        :param v2:
        :return:
        """
        barx1 = self.barx(v1)
        barx2 = self.barx(v2)
        if barx1 == barx2:
            return
        if barx1 >= len(self.histo) or barx2 >= len(self.histo):
            return
        self.histo[barx1] -= 1
        self.histo[barx2] += 1



if __name__ == "__main__":
    anim = GasAnimation()
    anim.start()
