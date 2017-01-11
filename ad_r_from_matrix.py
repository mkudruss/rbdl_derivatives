"""
Brief example on AD given in:

  Walter, S. F.:

    Structured higher-order algorithmic differentiation in the forward and
    reverse mode with application in optimum experimental design
"""

import os
import sys
import numpy as np

def roty(theta):
    return np.array([
        [np.cos(theta), 0.0, -np.sin(theta)],
        [      0.0,     1.0,            0.0],
        [np.sin(theta), 0.0,  np.cos(theta)],
    ])


def crm_3d(v):
    assert len(v) == 3 or v.shape == (3,)
    eye = np.eye(3)
    ret = np.zeros((3,3))
    ret[:,0] = np.cross(v, eye[0,:])
    ret[:,1] = np.cross(v, eye[1,:])
    ret[:,2] = np.cross(v, eye[2,:])
    return ret


def r_from_matrix(r, X):
    # cast as np.ndarrays
    r = np.asarray(r)
    X = np.asarray(X)

    # get E from X
    P_E = np.zeros((6, 3))
    P_E[:3, :] = np.eye(3)
    #print 'P_E\n', P_E
    E = P_E.transpose().dot(X).dot(P_E)
    #print 'E\n', E

    # get Erx from X
    P_Erx = np.zeros((6,3))
    P_Erx[3:, :] = np.eye(3)
    #print 'P_Erx\n', P_Erx
    Erx = -P_Erx.transpose().dot(X).dot(P_E)
    #print 'Erx\n', Erx

    # get r from Erx
    rx = E.transpose().dot(Erx)
    #print 'rx\n', rx
    eye = np.eye(3)
    P_rx0 = eye[:, 0:1]
    P_rx1 = eye[:, 1:2]
    P_rx2 = eye[:, 2:3]
    #print 'P_rx0:\n', P_rx0
    #print 'P_rx1:\n', P_rx1
    #print 'P_rx2:\n', P_rx2
    #print 'P_rx2.transpose().dot(rx.dot(P_rx1\n', P_rx2.transpose().dot(rx).dot(P_rx1)
    #print 'P_rx0.transpose().dot(rx.dot(P_rx2\n', P_rx0.transpose().dot(rx).dot(P_rx2)
    #print 'P_rx1.transpose().dot(rx.dot(P_rx0\n', P_rx1.transpose().dot(rx).dot(P_rx0)

    r0 = np.array([
        P_rx2.transpose().dot(rx).dot(P_rx1)[0, 0],
        P_rx0.transpose().dot(rx).dot(P_rx2)[0, 0],
        P_rx1.transpose().dot(rx).dot(P_rx0)[0, 0]
    ])
    print 'r0: ', r0

    # alternative r
    r1 = np.array([
        -P_rx1.transpose().dot(rx).dot(P_rx2)[0, 0],
        -P_rx2.transpose().dot(rx).dot(P_rx0)[0, 0],
        -P_rx0.transpose().dot(rx).dot(P_rx1)[0, 0]
    ])
    print 'r1: ', r1

    # assign outputs
    r[...] = r0


def r_from_matrix_dot(r, r_dot, X, X_dot):
    # cast as np.ndarrays
    r_dot = np.asarray(r_dot)
    X_dot = np.asarray(X_dot)
    r     = np.asarray(r)
    X     = np.asarray(X)
    P     = X_dot.shape[0]

    # get E from X
    P_E = np.zeros((6, 3))
    P_E[:3, :] = np.eye(3)
    #print 'P_E\n', P_E
    E_dot = [P_E.transpose().dot(x_dir).dot(P_E) for x_dir in X_dot[:]]
    E     = P_E.transpose().dot(X).dot(P_E)
    #print 'E\n', E

    # get Erx from X
    P_Erx = np.zeros((6,3))
    P_Erx[3:, :] = np.eye(3)
    #print 'P_Erx\n', P_Erx
    Erx_dot = [-P_Erx.transpose().dot(x_dir).dot(P_E) for x_dir in X_dot[:]]
    Erx     = -P_Erx.transpose().dot(X).dot(P_E)
    #print 'Erx\n', Erx

    # get r from Erx
    rx_dot = []
    for (E_dot_i, Erx_dot_i) in zip(E_dot, Erx_dot):
        rx_dot.append(E_dot_i.transpose().dot(Erx) + E.transpose().dot(Erx_dot_i))
    rx     = E.transpose().dot(Erx)
    #print 'rx\n', rx
    eye = np.eye(3)
    P_rx0 = eye[:, 0:1]
    P_rx1 = eye[:, 1:2]
    P_rx2 = eye[:, 2:3]
    #print 'P_rx0:\n', P_rx0
    #print 'P_rx1:\n', P_rx1
    #print 'P_rx2:\n', P_rx2
    #print 'P_rx2.transpose().dot(rx.dot(P_rx1\n', P_rx2.transpose().dot(rx).dot(P_rx1)
    #print 'P_rx0.transpose().dot(rx.dot(P_rx2\n', P_rx0.transpose().dot(rx).dot(P_rx2)
    #print 'P_rx1.transpose().dot(rx.dot(P_rx0\n', P_rx1.transpose().dot(rx).dot(P_rx0)

    r0_dot = [
        np.array([
            P_rx2.transpose().dot(rx_dot).dot(P_rx1)[0, 0],
            P_rx0.transpose().dot(rx_dot).dot(P_rx2)[0, 0],
            P_rx1.transpose().dot(rx_dot).dot(P_rx0)[0, 0]
        ]) for rx_dot in rx_dot[:]
    ]

    r0 = np.array([
        P_rx2.transpose().dot(rx).dot(P_rx1)[0, 0],
        P_rx0.transpose().dot(rx).dot(P_rx2)[0, 0],
        P_rx1.transpose().dot(rx).dot(P_rx0)[0, 0]
    ])
    #print 'r0: ', r0

    # alternative r
    r1 = np.array([
        -P_rx1.transpose().dot(rx).dot(P_rx2)[0, 0],
        -P_rx2.transpose().dot(rx).dot(P_rx0)[0, 0],
        -P_rx0.transpose().dot(rx).dot(P_rx1)[0, 0]
    ])
    #print 'r1: ', r1

    # assign outputs
    for i in range(P):
        r_dot[i,:] = r0_dot[i]
    r[...]     = r0


def fd(f, y, y_dot, x, x_dot, eps=1e-8):
    # cast as np.ndarrays
    y_dot = np.asarray(y_dot)
    x_dot = np.asarray(x_dot)
    y     = np.asarray(y)
    x     = np.asarray(x)

    # evaluation of function
    f(y, x)
    y_tmp = y.copy()  # save intermediate value

    # evaluation of directional derivative using finite differences
    P = x_dot.shape[0]
    for i in xrange(P):
        f(y_tmp, x + eps*x_dot[i, :])
        y_dot[i] = (y_tmp - y) / eps


if __name__ == '__main__':
    # set numpy printoptions
    np.set_printoptions(
        precision=18, threshold=None, edgeitems=None, linewidth=200,
        suppress=False, nanstr=None, infstr=None,
        #formatter={'float': lambda x: format(x, '1.08E')}
    )

    # setup input and output variables
    y = np.zeros((3,))

    # x = [ E    0 ]
    #     [-Erx  E ]
    x  = np.zeros((6,6))
    E  = roty(np.pi/3.)
    rx = crm_3d(np.array(range(1, 4)))
    print 'E\n',  E
    print 'rx\n', rx
    x[:3, :3] = E
    x[3:, :3] = -E.dot(rx)
    x[3:, 3:] = E
    print 'x:\n', x

    # setup derivatives and directions
    P = np.prod(x.shape)  # number of directions
    y_dot = np.zeros(((P,) + y.shape))  # directional derivative (y.shape, P)
    # NOTE: passing the identity will result in computing the Jacobian
    x_dot = np.zeros(((P,) + x.shape) )    # directions (x.shape, P)
    for i in range(6):
        for j in range(6):
            x_dot[i*6 + j, i, j] = 1.0
    print x_dot.shape

    # function evaluation
    print 'FUNCTION EVALUATION'
    print 'x = \n', x

    r_from_matrix(y, x)
    print 'y = f(x) =\n', y

    # derivative evaluation
    print 'DERIVATIVE EVALUATION'
    print 'x_dot = \n', x_dot

    r_from_matrix_dot(y, y_dot, x, x_dot)
    # save reference values
    y_r = y.copy()
    y_dot_r = y_dot.copy()
    print 'y_dot = f_dot(x) =\n', y_dot.shape

    fd(r_from_matrix, y, y_dot, x, x_dot, eps=1e-08)
    print 'y_dot = f_fd(x) =\n', y_dot
    print 'err =\n', np.abs(y_dot - y_dot_r)
