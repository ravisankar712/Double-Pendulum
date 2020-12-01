''' Author : Ravisankar R -- https://github.com/ravisankar712 '''
''' Inspired from -- https://www.youtube.com/watch?v=d0Z8wLLPNE0&t=24s '''

from manimlib.imports import *

#gravity
g = 9.81

class Pendulum(VGroup):
    CONFIG = {
        "m" : 1.0,
        "L" : 1.0,
        "start_angle" : 0.0,
        "pivot" : ORIGIN + 2 * UP,
        "bob_color" : BLUE
    }
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
        self.angle = self.start_angle
        self.omega = 0.0
        self.create_body()

    def trigger_motion(self):
        self.add_updater(lambda m, dt: m.oscillation(dt))
    
    def create_body(self):
        #assumes origin as the pivot
        end = self.L * np.array([np.cos(self.angle - PI/2), np.sin(self.angle - PI/2), 0.0])
        rod = Line(self.pivot, self.pivot + end, stroke_width=8.0)
        rod.set_color(GREY)
        rod.move_to(self.pivot + end/2.0)
        bob = Circle(radius=0.2, fill_opacity=1.0).set_color(self.bob_color)
        bob.move_to(self.pivot + end)
        self.rod = rod
        self.bob = bob
        self.add(self.rod, self.bob)

    def oscillation(self, dt):
        #Euler Integration
        acc = - (g / self.L) * np.sin(self.angle)
        self.omega += acc * dt
        self.angle += self.omega * dt
        end = self.L * np.array([np.cos(self.angle - PI/2), np.sin(self.angle - PI/2), 0.0])
        self.bob.move_to(self.pivot + end)
        self.rod.put_start_and_end_on(self.pivot, self.pivot + end)

class Path(VGroup):

    CONFIG = {
        "trail_length" : 20
    }

    def __init__(self, bob, col=WHITE, **kwargs):
        super().__init__(**kwargs)
        self.col = col
        self.bob = bob
        # self.offset = pendulum.pivot
        self.points = [bob.get_center()]
        self.path = VGroup()
        self.add(self.path)
        self.add_updater(lambda m, dt: m.draw_path(dt))

    def draw_path(self, dt):
        l = len(self.points)
        if l > self.trail_length:
            self.points.pop(0)
        pt = self.bob.get_center()
        if np.linalg.norm(pt - self.points[-1]) != 0:
            self.points.append(pt)
        current_path = VGroup()
        for i in range(l - 1):
            line = Line(self.points[i], self.points[i+1], color=self.col)
            current_path.add(line)
        self.path.become(current_path)

class DoublePendulum(VGroup):
    CONFIG = {
        "m1" : 1.0,
        "m2" : 1.0,
        "L1" : 2.0,
        "L2" : 2.0,
        "pivot" : ORIGIN + UP * 2,
        "bob1_color" : BLUE,
        "bob2_color" : GREEN
    }

    def __init__(self, theta1=0.0, theta2=0.0, **kwargs):
        super().__init__(**kwargs)
        self.theta1 = theta1
        self.theta2 = theta2
        self.omega1 = 0.0
        self.omega2 = 0.0

        self.create_body()

    def trigger_motion(self):
        self.add_updater(lambda m, dt: m.oscillation(dt))

    def create_body(self):
        #assumes origin as the pivot
        end1 = self.L1 * np.array([np.cos(self.theta1 - PI/2), np.sin(self.theta1 - PI/2), 0.0])
        rod1 = Line(self.pivot, self.pivot + end1, stroke_width=8.0)
        rod1.set_color(GREY)
        rod1.move_to(self.pivot + end1/2.0)
        bob1 = Circle(radius=0.2, fill_opacity=1.0).set_color(self.bob1_color)
        bob1.move_to(self.pivot + end1)

        end2 = self.L2 * np.array([np.cos(self.theta2 - PI/2), np.sin(self.theta2 - PI/2), 0.0])
        rod2 = Line(self.pivot + end1, self.pivot + end1 + end2, stroke_width=8.0)
        rod2.set_color(GREY)
        rod2.move_to(self.pivot + end1 + end2/2.0)
        bob2 = Circle(radius=0.2, fill_opacity=1.0).set_color(self.bob2_color)
        bob2.move_to(self.pivot + end1 + end2)
        self.rod1 = rod1
        self.bob1 = bob1
        self.rod2 = rod2
        self.bob2 = bob2
        self.add(self.rod1, self.rod2, self.bob2, self.bob1)

    def oscillation(self, dt):
        m1 = self.m1
        m2 = self.m2
        L1 = self.L1
        L2 = self.L2
        theta1 = self.theta1
        theta2 = self.theta2
        omega1 = self.omega1
        omega2 = self.omega2
        #The Euler Method. Update to Runge Kutta for more accuracy.
        den = L1 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2))
        acc1 = (-g * (2 * m1 + m2) * np.sin(theta1) + 
                -m2 * g * np.sin(theta1 - 2 * theta2) + 
                -2 * np.sin(theta1 - theta2) * m2 * omega2**2 * L2 + omega1**2 * L1 * np.cos(theta1 - theta2)) / den

        den = L2 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2))
        acc2 = (2 * np.sin(theta1 - theta2) * 
                (omega1**2 * L1 * (m1 + m2) + g * (m1 + m2) * np.cos(theta1) + omega2**2 * L2 * m2 * np.cos(theta1 - theta2))) / den
        
        self.omega1 += acc1 * dt
        self.theta1 += self.omega1 * dt
        end1 = self.L1 * np.array([np.cos(self.theta1 - PI/2), np.sin(self.theta1 - PI/2), 0.0])
        self.bob1.move_to(self.pivot + end1)
        self.rod1.put_start_and_end_on(self.pivot, self.pivot + end1)

        self.omega2 += acc2 * dt
        self.theta2 += self.omega2 * dt
        end2 = self.L2 * np.array([np.cos(self.theta2 - PI/2), np.sin(self.theta2 - PI/2), 0.0])
        self.bob2.move_to(self.pivot + end1 + end2)
        self.rod2.put_start_and_end_on(self.pivot + end1, self.pivot + end1 + end2)
        

class Test(Scene):
    def construct(self):
        # p = Pendulum(L=2.0, start_angle=PI/2 + PI/4, bob_color=GREEN)
        # path = Path(p.bob, p.bob_color)
        # self.add(p, path)
        dp = DoublePendulum(theta1 = PI/2, theta2=PI/2, bob1_color="#ff1654", bob2_color="#247ba0")
        p1 = Path(dp.bob1, dp.bob1_color)
        p2 = Path(dp.bob2, dp.bob2_color)
        self.add(dp, p1, p2)
        self.wait(30)

class Pend_Intro(Scene):
    def construct(self):
        p = Pendulum(L = 2.0, bob_color="#f94144")
        self.add(p)
        self.wait()
        p.generate_target()
        p.target = Pendulum(L = 2.0, bob_color="#f94144", start_angle=PI/3)

        self.play(
            MoveToTarget(p)
        )

        self.wait()
        self.remove(p)
        p = Pendulum(L = 2.0, bob_color="#f94144", start_angle=PI/3)
        p.trigger_motion()
        path = Path(p.bob, p.bob_color)
        self.add(p, path)
        self.wait(12)

class Pend_with_different_starts(Scene):
    def construct(self):
        pends = []
        colors = ["#f94144", "#f3722c", "#f8961e", "#f9c74f", "#90be6d", "#43aa8b", "#577590"]
        ang = PI/3
        dtheta = 0.1/len(colors)
        for i in range(len(colors)):
            p = Pendulum(L = 2.0, bob_color=colors[i], start_angle=ang+i*dtheta)
            pends.append(p)

        # self.add(*pends)
        # self.wait()
        self.play(
            AnimationGroup(
                *[FadeIn(p) for p in pends], lag_ratio=0.5
            )
        )
        self.wait()
        for p in pends:
            self.add(Path(p.bob, p.bob_color))
            p.trigger_motion()
        self.wait(20)

class DBPend_Intro(Scene):
    def construct(self):
        p = DoublePendulum(L1 = 2.0, bob1_color="#f94144", L2 = 2.0, bob2_color="#90be6d")
        self.add(p)
        self.wait()
        p.generate_target()
        p.target =  DoublePendulum(theta1 = PI/2, theta2=PI/2, L1 = 2.0, bob1_color="#f94144", L2 = 2.0, bob2_color="#90be6d")

        self.play(
            MoveToTarget(p)
        )

        self.wait()
        self.remove(p)
        p = DoublePendulum(theta1 = PI/2, theta2=PI/2, L1 = 2.0, bob1_color="#f94144", L2 = 2.0, bob2_color="#90be6d")
        p.trigger_motion()
        path1 = Path(p.bob1, p.bob1_color)
        path2 = Path(p.bob2, p.bob2_color)
        self.add(p, path1, path2)
        self.wait(30)

class DBPend_with_different_starts(Scene):
    def construct(self):
        pends = []
        colors = ["#f94144", "#f3722c", "#f8961e", "#f9c74f", "#90be6d", "#43aa8b", "#577590"]
        # colors = ["#023e8a", "#00b4d8", "#d00000", "#dc2f02", "#90be6d", "#43aa8b", "#f0ead2", "#dde5b6", "#3c096c", "#7b2cbf", "#ff97b7", "#ff5d8f"]
        ang = PI/2
        dtheta = 0.1/len(colors)
        for i in range(len(colors)):
            p = DoublePendulum(theta1 = PI/2, theta2=ang + i*dtheta, L1 = 2.0, bob1_color=colors[i], L2 = 2.0, bob2_color=colors[i])
            pends.append(p)

        # self.add(*pends)
        # self.wait()
        self.play(
            AnimationGroup(
                *[FadeIn(p) for p in pends], lag_ratio=0.5
            )
        )
        self.wait()

        for p in pends:
            self.add(Path(p.bob2, p.bob2_color, trail_length=50))
            p.trigger_motion()
        self.wait(35)


class Intro(Scene):
    def construct(self):
        sciencesort = Text('science.sort', font='Lucida Console')
        # dot = TexMobject("\\cdot")
        # sort = Text('sort', font='Lucida Console')
        langle = TexMobject("\\langle").set_color(RED)
        rangle = TexMobject("\\rangle").set_color(RED)
        
        sciencesort.next_to(langle, RIGHT, buff=0.1)
        # dot.next_to(science, RIGHT, buff=0.0)
        # sort.next_to(dot, RIGHT, buff=0.0)
        rangle.next_to(sciencesort, RIGHT, buff=0.1)

        logo = VGroup(langle, sciencesort, rangle).move_to(ORIGIN)
        self.play(
            FadeInFromDown(logo)
        )
        self.wait(0.5)
        langle1 = TexMobject("\\langle").set_color(RED)
        rangle1 = TexMobject("\\rangle").set_color(RED)
        scienceshort = Text('science.short', font='Lucida Console').next_to(langle1, RIGHT, buff=0.1)
        rangle1.next_to(scienceshort, RIGHT, buff=0.1)
        logo1 = VGroup(langle1, scienceshort, rangle1).move_to(ORIGIN)
        
        self.play(
            Transform(logo, logo1)
        )

        self.wait()

SCENES_IN_ORDER = [
Intro, 
Pend_Intro, 
Pend_with_different_starts, 
DBPend_Intro, 
DBPend_with_different_starts
]