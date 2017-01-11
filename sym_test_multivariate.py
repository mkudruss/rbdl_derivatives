# enable pretty printing
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

def Xm_mat_E_mErx (E, mErx):
    return Matrix ((
        (E[0,0], E[0,1], E[0,2], 0, 0, 0,),
        (E[1,0], E[1,1], E[1,2], 0, 0, 0,),
        (E[2,0], E[2,1], E[2,2], 0, 0, 0,),
        (mErx[0,0], mErx[0,1], mErx[0,2], E[0,0], E[0,1], E[0,2]),
        (mErx[1,0], mErx[1,1], mErx[1,2], E[1,0], E[1,1], E[1,2]),
        (mErx[2,0], mErx[2,1], mErx[2,2], E[2,0], E[2,1], E[2,2])
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
    mErx = S_Erx.transpose() * X * S_E
    rx = -E.transpose() * mErx

    S_rx_0 = Matrix (( (1), (0), (0) ))
    S_rx_1 = Matrix (( (0), (1), (0) ))
    S_rx_2 = Matrix (( (0), (0), (1) ))

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

    # Expression of E with individual entries
    E_00 = Symbol('E(0,0)')
    E_01 = Symbol('E(0,1)')
    E_02 = Symbol('E(0,2)')
    E_10 = Symbol('E(1,0)')
    E_11 = Symbol('E(1,1)')
    E_12 = Symbol('E(1,2)')
    E_20 = Symbol('E(2,0)')
    E_21 = Symbol('E(2,1)')
    E_22 = Symbol('E(2,2)')

    E = Matrix ((
        (E_00, E_01, E_02),
        (E_10, E_11, E_12),
        (E_20, E_21, E_22),
        ))

    # Expression of r with individual values
    r_0 = Symbol('r[0]')
    r_1 = Symbol('r[1]')
    r_2 = Symbol('r[2]')

    r = Matrix (( r_0, r_1, r_2 ))

    # Expression of rx using individual values
    rx = Matrix ((
        (0, -r[2], r[1]),
        (r[2], 0, -r[0]),
        (-r[1], r[0], 0)
        ))

    Erx = E * rx

    # Expression of -Erx with individual values for the entries
    Erx_00 = Symbol('Erx(0,0)')
    Erx_01 = Symbol('Erx(0,1)')
    Erx_02 = Symbol('Erx(0,2)')

    Erx_10 = Symbol('Erx(1,0)')
    Erx_11 = Symbol('Erx(1,1)')
    Erx_12 = Symbol('Erx(1,2)')

    Erx_20 = Symbol('Erx(2,0)')
    Erx_21 = Symbol('Erx(2,1)')
    Erx_22 = Symbol('Erx(2,2)')

    mErx = -Matrix ((
        (Erx_00, Erx_01, Erx_02),
        (Erx_10, Erx_11, Erx_12),
        (Erx_20, Erx_21, Erx_22),
        ))

    print ("mErx")
    show (mErx)

    print ("- E * rx = ")
    show (- E * rx)
    manual_mErx = - E * rx
    print ("- E * rx - mErx")
    show ((- E * rx - mErx).subs ([
        # replace mErx matrix with symbolic values
        (Erx_00, Erx[0,0]), (Erx_01, Erx[0,1]), (Erx_02, Erx[0,2]),
        (Erx_10, Erx[1,0]), (Erx_11, Erx[1,1]), (Erx_12, Erx[1,2]),
        (Erx_20, Erx[2,0]), (Erx_21, Erx[2,1]), (Erx_22, Erx[2,2]),

        (r_0,1), (r_1,2), (r_2,3),

        (E_00,1), (E_01,0), (E_02,0),
        (E_10,0), (E_11,1), (E_12,0),
        (E_20,0), (E_21,0), (E_22,1),

        ]
        ))

    #    X = Xm (E, Matrix ((Symbol('r[0]'), Symbol('r[1]'), Symbol('r[2]'))) * q0)
    X = Xm_mat_E_mErx (E, mErx * q0)
    print ("d X / dq0")
    show (X.diff(q0))
    print ("d Xm_r /dq")
    show (Xm_calc_r(X).diff(q0))
    print (values)
    show (Xm_calc_r(X).diff(q0).subs([
        # replace mErx matrix with symbolic values
        (Erx_00, Erx[0,0]), (Erx_01, Erx[0,1]), (Erx_02, Erx[0,2]),
        (Erx_10, Erx[1,0]), (Erx_11, Erx[1,1]), (Erx_12, Erx[1,2]),
        (Erx_20, Erx[2,0]), (Erx_21, Erx[2,1]), (Erx_22, Erx[2,2]),

        (r_0,1), (r_1,2), (r_2,3),

        (E_00,1), (E_01,0), (E_02,0),
        (E_10,0), (E_11,1), (E_12,0),
        (E_20,0), (E_21,0), (E_22,1),

        ]))

#    print ("X:")
#    show (X)
#
#    print ("dX/dq0:")
#    show (X.diff(q0))
#
#    print ("dXm_calc_r (X)/dq0:")
#    r = Xm_calc_r (X)
#    show (r.diff (q0))


