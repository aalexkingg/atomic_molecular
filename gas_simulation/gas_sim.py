from vpython import *
import base.atomic as atomic

# constants
WIDTH = 800
HEIGHT = 400

k = 1.4E-23  # Boltzmann constant


class GasAnimation:
    def __init__(self):
        """

        """
        self.deltav = 100  # binning for v histogram
        self.nhisto = 0  # number of histogram snapshots to average

        self.animation = canvas(width=WIDTH / 3, height=HEIGHT, align='left')
        self.animation.userpan = True

        # self.particle_count_slider = slider(bind=set_count, min=1, max=300, align='right')

        self.draw_speed_graph()

        self.num_atoms = 20  # change this to have more or fewer atoms

        # Typical values
        self.L = 1  # container is a cube L on a side
        gray = color.gray(0.7)  # color of edges of container
        self.mass = 4E-3 / 6E23  # helium mass
        self.Ratom = 0.03  # wildly exaggerated size of helium atom

        self.T = 300  # around room temperature
        self.dt = 1E-5




    def start(self):
        """

        :return:
        """
        Atoms = []
        p = []
        apos = []
        pavg = sqrt(2 * mass * 1.5 * k * T)  # average kinetic energy p**2/(2mass) = (3/2)kT

        nhisto = int(4500 / deltav)
        histo = []
        for i in range(nhisto): histo.append(0.0)
        histo[barx(pavg / mass)] = Natoms

        accum = []
        for i in range(int(3000 / deltav)): accum.append([deltav * (i + .5), 0])
        vdist = gvbars(color=color.red, delta=deltav)

        while True:
            rate(300)
            # Accumulate and average histogram snapshots
            for i in range(len(accum)): accum[i][1] = (nhisto * accum[i][1] + histo[i]) / (nhisto + 1)
            if nhisto % 10 == 0:
                vdist.data = accum
            nhisto += 1

            # Update all positions
            for i in range(Natoms): Atoms[i].pos = apos[i] = apos[i] + (p[i] / mass) * dt

            # Check for collisions
            hitlist = self.checkCollisions()

            # If any collisions took place, update momenta of the two atoms
            for ij in hitlist:
                i = ij[0]
                j = ij[1]
                ptot = p[i] + p[j]
                posi = apos[i]
                posj = apos[j]
                vi = p[i] / mass
                vj = p[j] / mass
                vrel = vj - vi
                a = vrel.mag2
                if a == 0: continue;  # exactly same velocities
                rrel = posi - posj
                if rrel.mag > Ratom: continue  # one atom went all the way through another

                # theta is the angle between vrel and rrel:
                dx = dot(rrel, vrel.hat)  # rrel.mag*cos(theta)
                dy = cross(rrel, vrel.hat).mag  # rrel.mag*sin(theta)
                # alpha is the angle of the triangle composed of rrel, path of atom j, and a line
                #   from the center of atom i to the center of atom j where atome j hits atom i:
                alpha = asin(dy / (2 * Ratom))
                d = (2 * Ratom) * cos(alpha) - dx  # distance traveled into the atom from first contact
                deltat = d / vrel.mag  # time spent moving from first contact to position inside atom

                posi = posi - vi * deltat  # back up to contact configuration
                posj = posj - vj * deltat
                mtot = 2 * mass
                pcmi = p[i] - ptot * mass / mtot  # transform momenta to cm frame
                pcmj = p[j] - ptot * mass / mtot
                rrel = norm(rrel)
                pcmi = pcmi - 2 * pcmi.dot(rrel) * rrel  # bounce in cm frame
                pcmj = pcmj - 2 * pcmj.dot(rrel) * rrel
                p[i] = pcmi + ptot * mass / mtot  # transform momenta back to lab frame
                p[j] = pcmj + ptot * mass / mtot
                apos[i] = posi + (p[i] / mass) * deltat  # move forward deltat in time
                apos[j] = posj + (p[j] / mass) * deltat
                interchange(vi.mag, p[i].mag / mass)
                interchange(vj.mag, p[j].mag / mass)

            for i in range(Natoms):
                loc = apos[i]
                if abs(loc.x) > L / 2:
                    if loc.x < 0:
                        p[i].x = abs(p[i].x)
                    else:
                        p[i].x = -abs(p[i].x)

                if abs(loc.y) > L / 2:
                    if loc.y < 0:
                        p[i].y = abs(p[i].y)
                    else:
                        p[i].y = -abs(p[i].y)

                if abs(loc.z) > L / 2:
                    if loc.z < 0:
                        p[i].z = abs(p[i].z)
                    else:
                        p[i].z = -abs(p[i].z)

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
        d = L / 2 + Ratom
        r = 0.005
        boxbottom = curve(color=gray, radius=r)
        boxbottom.append([vector(-d, -d, -d), vector(-d, -d, d), vector(d, -d, d), vector(d, -d, -d), vector(-d, -d, -d)])
        boxtop = curve(color=gray, radius=r)
        boxtop.append([vector(-d, d, -d), vector(-d, d, d), vector(d, d, d), vector(d, d, -d), vector(-d, d, -d)])
        vert1 = curve(color=gray, radius=r)
        vert2 = curve(color=gray, radius=r)
        vert3 = curve(color=gray, radius=r)
        vert4 = curve(color=gray, radius=r)
        vert1.append([vector(-d, -d, -d), vector(-d, d, -d)])
        vert2.append([vector(-d, -d, d), vector(-d, d, d)])
        vert3.append([vector(d, -d, d), vector(d, d, d)])
        vert4.append([vector(d, -d, -d), vector(d, d, -d)])

    def create_particles(self):
        """

        :return:
        """
        for i in range(Natoms):
            x = L * random() - L / 2
            y = L * random() - L / 2
            z = L * random() - L / 2
            if i == 0:
                Atoms.append(
                    sphere(pos=vector(x, y, z), radius=Ratom / 2, color=color.cyan, make_trail=False, retain=100,
                           trail_radius=0.3 * Ratom))
            else:
                Atoms.append(sphere(pos=vector(x, y, z), radius=Ratom / 2, color=gray))
            apos.append(vec(x, y, z))
            theta = pi * random()
            phi = 2 * pi * random()
            px = pavg * sin(theta) * cos(phi)
            py = pavg * sin(theta) * sin(phi)
            pz = pavg * cos(theta)
            p.append(vector(px, py, pz))

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
            ymax=1000
        )
        theory = gcurve(color=color.blue, width=2)

        dv = 10
        for v in range(0, 3001 + dv, dv):  # theoretical prediction
            theory.plot(v, (deltav / dv) * Natoms * 4 * pi * ((mass / (2 * pi * k * T)) ** 1.5) * exp(
                -0.5 * mass * (v ** 2) / (k * T)) * (v ** 2) * dv)

    def checkCollisions(self):
        """

        :return:
        """
        hitlist = []
        r2 = 2 * Ratom
        r2 *= r2
        for i in range(Natoms):
            ai = apos[i]
            for j in range(i):
                aj = apos[j]
                dr = ai - aj
                if mag2(dr) < r2: hitlist.append([i, j])
        return hitlist




    def barx(self, v):
        """

        :param v:
        :return:
        """
        return int(v / deltav)  # index into bars array


    def set_count(self, evt):
        """

        :param evt:
        :return:
        """
        global Natoms
        Natoms = floor(evt.value)


    def interchange(self, v1, v2):  # remove from v1 bar, add to v2 bar
        """

        :param v1:
        :param v2:
        :return:
        """
        barx1 = barx(v1)
        barx2 = barx(v2)
        if barx1 == barx2:  return
        if barx1 >= len(histo) or barx2 >= len(histo): return
        histo[barx1] -= 1
        histo[barx2] += 1



if __name__ == "__main__":
    anim = GasAnimation()
    anim.start()
