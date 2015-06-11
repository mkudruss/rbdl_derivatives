 # ennable pretty printing
from sympy import pprint

# sympy linear algebra
from sympy.matrices import *
from sympy import Symbol, Eq, Expr
from sympy import simplify, nsimplify, trigsimp, expand
from sympy import collect
from sympy import diff
from sympy import solve
from sympy import S, pi
from sympy import sin, cos


def _ry(theta=Symbol('theta')):
    s = sin(theta)
    c = cos(theta)
    return Matrix((
        (   c,    0,   -s),
        (   0,    1,    0),
        (   s,    0,    c),
    ))


def crm_3d(v):
    assert len(v) == 3 or v.shape == (3,1)
    v = Matrix(len(v), 1, v)
    ret = zeros(3,3)
    ret[:,0] = v.cross(Matrix(3,1,[1, 0, 0]))
    ret[:,1] = v.cross(Matrix(3,1,[0, 1, 0]))
    ret[:,2] = v.cross(Matrix(3,1,[0, 0, 1]))
    return ret

def Xm (E=Symbol('E'), r=Symbol('r')):
    Erx = -E * crm_3d (r)
    return Matrix ((
        (E[0,0], E[0,1], E[0,2], 0, 0, 0,),
        (E[1,0], E[1,1], E[1,2], 0, 0, 0,),
        (E[2,0], E[2,1], E[2,2], 0, 0, 0,),
        (-Erx[0,0], -Erx[0,1], -Erx[0,2], E[0,0], E[0,1], E[0,2]),
        (-Erx[1,0], -Erx[1,1], -Erx[1,2], E[1,0], E[1,1], E[1,2]),
        (-Erx[2,0], -Erx[2,1], -Erx[2,2], E[2,0], E[2,1], E[2,2])
        ))

def Xtrans (r=Matrix ((Symbol('r[0]'), Symbol('r[1]'), Symbol('r[2]')))):
    return Xm (eye(3), r)

def Xm_calc_r (X):
    S_E = Matrix ((
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        (0, 0, 0),
        (0, 0, 0),
        (0, 0, 0),
        ))

    S_Erx = Matrix ((
        (0, 0, 0),
        (0, 0, 0),
        (0, 0, 0),
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        ))

    E = S_E.transpose() * X * S_E
    Erx = S_Erx.transpose() * X * S_E
    rx = E.transpose() * Erx

    S_rx_0 = Matrix ((
        (1), (0), (0)
        ))
    S_rx_1 = Matrix ((
        (0), (1), (0)
        ))
    S_rx_2 = Matrix ((
        (0), (0), (1)
        ))


    r = Matrix ((
        S_rx_2.transpose() * rx * S_rx_1,
        S_rx_0.transpose() * rx * S_rx_2,
        S_rx_1.transpose() * rx * S_rx_0
        ))

    # alternative r
    r_2 = Matrix ((
        -S_rx_1.transpose() * rx * S_rx_2,
        -S_rx_2.transpose() * rx * S_rx_0,
        -S_rx_0.transpose() * rx * S_rx_1
        ))

    return r

def show(expr):
    print pprint(expr)

if __name__ == '__main__':
    q0 = Symbol('q_0')
    q1 = Symbol('q_1')

    sdict = {q0 : pi/2, q1 : pi/3}

    A = _ry(q0)
    B = _ry(q1)

    y = A*B
    show(y)
    show(y.subs(sdict))

    y_q0 = y.diff(q0)
    y_q1 = y.diff(q1)
    show(y_q0)
    show(y_q0.subs(sdict))
    show(y_q1)
    show(y_q1.subs(sdict))

    X = Xtrans (Matrix ((1, 0, 0)) * q0)
    show (X)
    show (X.diff(q0))

    r = Xm_calc_r (X)
    show (r.diff (q0))
    
