# standard preamble
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
from Timba import pynode
from Timba import dmi
from Timba import utils
from Timba.Meq import meq

from numpy import *
from numpy import linalg


_dbg = utils.verbosity(0, name='test_pynode')
_dprint = _dbg.dprint
_dprintf = _dbg.dprintf

# This class is meant to illustrate the pynode interface. All pynodes
# need to be derived from the Timba.pynode.PyNode class.


class KLNode (pynode.PyNode):

    '''creates the U matrix on piercepoints defined by the children, with rank = rank'''
    # An __init__ method is not necessary at all. If you do define your
    # own for some reason, make sure you call pynode.PyNode.__init__()
    # as done below

    def __init__(self, *args):
        pynode.PyNode.__init__(self, *args)
        self.set_symdeps('domain')

        # Finally, we should define a get_result method. (We don't have to, but if
        # we don't, then what's the point of this node?)
        # This is called with a list of child results (possibly empty, if we're a
        # leaf node), and a request object
    def get_result(self, request, *children):
        if len(children) < 2:
            raise TypeError("this is NOT a leaf node, At least 2  children expected!")
        cells = request.cells
        res = meq.result(None, cells)

        Xp_table = [0, 0, 0]
        # if necessary loop over time and other axis
        # First get PiercePoints
        res1 = children[1]
        vs1 = res1.vellsets
        # num_pierce_points each containing vector of length 2 or 3 (x,y(,z))
        dims = res1.dims
        vector_size = dims[1]
        num_pp = dims[0]
        x = []
        for j in range(num_pp):
            x.append([])
            for i in range(vector_size):
                x[j].append(0)

        # x=[[0,]*vector_size,]*num_pp;
        Xp_table = [0, ] * num_pp
        for i in range(num_pp):
            for j in range(vector_size):
                x[i][j] = vs1[i * vector_size + j].value
        shape = meq.shape(x[0][0])
        grid_size = len(x[0][0])
        # print "grid_size",grid_size;

        # Then the parameters
        res0 = children[0]
        vs0 = res0.vellsets
        # should be one per parameter
        num_parms = len(vs0)
        val0 = [0, ] * num_parms
        par = [0, ] * num_parms
        solvable = False
        parmidx = list(range(num_parms))

        # initialize  result vector 1 field per pp;
        #      ph=[[0,]*grid_size,]*num_parms;
        ph = []
        for j in range(num_pp):
            ph.append([])
            for i in range(grid_size):
                ph[j].append(0)

        # if solvable they should have perturbed_values;
        if 'perturbed_value' in vs0[0]:
            solvable = True
            pert = [0] * num_parms
            pt = []
            perturbations = ()
            spid_index = ()
            # get spidindex
            for i in range(num_parms):
                spid_index += vs0[i].spid_index

            # sort on spidindex
            for i in range(num_parms):
                for j in range(i):
                    if spid_index[i] < spid_index[parmidx[j]]:
                        newidx = i
                        for k in range(j, i + 1):
                            old = parmidx[k]
                            parmidx[k] = newidx
                            newidx = old
                        break
            spid_index = ()
            # 1 perturbed value per pp per paramter
            # or the other way around???
            for i in range(num_parms):
                pt.append([])

                for i2 in range(num_pp):
                    pt[i].append([])
                    for g in range(grid_size):
                        pt[i][i2].append(0)
        # fill values sorted
        for pn in range(num_parms):
            pnidx = parmidx[pn]
            val0[pn] = vs0[pnidx].value
            par[pn] = val0[pn][0]
            # use first value (in case of constant parms)

            # assume single perturbation parms for the moment
            if solvable:
                pert[pn] = vs0[pnidx].perturbed_value[0]
                spid_index += vs0[pnidx].spid_index
                perturbations += vs0[pnidx].perturbations

        # get U matrix
        # Only if domain changed!!
        for i in range(grid_size):
            for j in range(vector_size):
                for k in range(num_pp):
                    Xp_table[k] = x[k][j][i]
            U, S, Ut, B = get_U(Xp_table, order=num_parms)

            for pn in range(num_parms):
                if len(val0[pn]) > i:
                    # incorrect check for freq/time dependendies
                    par[pn] = (val0[pn][i])

            value = dot(U, par)
            for v in range(num_pp):
                ph[v][i] = value[v]
            for v in range(num_parms):
                if solvable:
                    old = par[v]
                    if len(pert[v]) > i:
                        par[v] = pert[v][i]
                    else:
                        par[v] = pert[v][0]
                    pval = dot(U, par)
                    par[v] = old
                    for v2 in range(num_pp):
                        pt[v][v2][i] = pval[v2]

        vellsets = []
        for v in range(num_pp):
            value = meq.vells(shape)
            ph[v] = array(ph[v])
            ph[v].shape = shape
            value[:] = ph[v]

            vellsets.append(meq.vellset(value))
            if solvable:
                vellsets[v].perturbed_value = ()
                for v2 in range(num_parms):
                    value = meq.vells(shape)
                    pt[v2][v] = array(pt[v2][v])
                    pt[v2][v].shape = shape
                    value[:] = pt[v2][v]
                    vellsets[v].perturbed_value += (value,)
                    vellsets[v].spid_index = spid_index
                    vellsets[v].perturbations = perturbations
        res.vellsets = vellsets
        # res.dims = [3,1];
        #    return meq.result(meq.vellset(value),cells);

        # print "result",res
        return res


def get_U(Xp_table=None, order=10, beta=3. / 5., r_0=1):
        # bepaling van U

        # calculate structure matrix
    p_count = len(Xp_table)
    D_table = resize(Xp_table, (p_count, p_count, 2))
    D_table = transpose(D_table, (1, 0, 2)) - D_table
    D_table = add.reduce(D_table ** 2, 2)
    D_table = (D_table / (r_0 ** 2)) ** (beta / 2.)

    # calculate covariance matrix C
    # calculate partial product for interpolation B
    C_table = - D_table / 2.
    C_table = transpose(C_table - (add.reduce(C_table, 0) / float(p_count)))
    B_table = add.reduce(C_table, 0) / float(p_count)
    C_table = C_table - B_table

    # eigenvalue decomposition
    # select subset of base vectors
    [U_table, S, Ut_table] = linalg.svd(C_table)
    U_table = U_table[:, 0: order]
    S = S[0: order]
    Ut_table = Ut_table[0: order, :]
    Si = 1. / S

    return [U_table, Si, Ut_table, B_table]


def get_interpol(Xp_table, U_table, P, Si, Ut_table):
# interpolatie naar nieuwe X

        # extract gradient
    Q = dot(transpose(Xp_table), Xp_table)
    Q = dot(linalg.inv(Q), transpose(Xp_table))
    Q = dot(Q, U_table)
    Q = dot(Q, transpose([P]))
    Q = (transpose(Q))[0]
    R = dot([Q], transpose(Xp_table))
    R = (dot(R, U_table))[0]
    R = P - R

    # calculate interpolation matrix
    F_table = transpose(dot([Si * R], Ut_table))

    PHI = phi_klmap_model(X_table, Xp_table, B_table, F_table, beta=beta, r_0=r_0)


def phi_klmap_model(X, Xp_table, B_table, F_table, beta=5. / 3., r_0=1.):
# B_table = (1/m)(1T)(C_table)(A_table)
# F_table = ( Ci_table )( U_table )( P_table )

    # input check
    if (len(shape(X)) == 1):
        X_table = array([X])
    elif (len(shape(X)) == 2):
        X_table = X

    # calculate structure matrix
    x_count = len(X_table)
    p_count = len(Xp_table)
    D_table = transpose(resize(X_table, (p_count, x_count, 2)), (1, 0, 2))
    D_table = D_table - resize(Xp_table, (x_count, p_count, 2))
    D_table = add.reduce(D_table ** 2, 2)
    D_table = (D_table / (r_0 ** 2)) ** (beta / 2.)

    # calculate covariance matrix
    C_table = - D_table / 2.
    C_table = transpose(transpose(C_table) - (add.reduce(C_table, 1) / float(p_count)))
    C_table = C_table - B_table

    phi = dot(C_table, F_table)
    phi = reshape(phi, (x_count))
    if (len(phi) == 1):
        phi = phi[0]

    return phi
