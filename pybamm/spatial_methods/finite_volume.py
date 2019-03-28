#
# Finite Volume discretisation class
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm

import numpy as np
from scipy.sparse import spdiags
from scipy.sparse import eye
from scipy.sparse import kron


class FiniteVolume(pybamm.SpatialMethod):
    """
    A class which implements the steps specific to the finite volume method during
    discretisation.

    Parameters
    ----------
    mesh : :class:`pybamm.Mesh` (or subclass)
        Contains all the submeshes for discretisation

    **Extends:"": :class:`pybamm.SpatialMethod`
    """

    def __init__(self, mesh):
        # add npts_for_broadcast to mesh domains for this particular discretisation
        for dom in mesh.keys():
            for i in range(len(mesh[dom])):
                mesh[dom][i].npts_for_broadcast = mesh[dom][i].npts
        super().__init__(mesh)

    def spatial_variable(self, symbol):
        """
        Creates a discretised spatial variable compatible with
        the FiniteVolume method.

        Parameters
        -----------
        symbol : :class:`pybamm.SpatialVariable`
            The spatial variable to be discretised.

        Returns
        -------
        :class:`pybamm.Vector`
            Contains the discretised spatial variable
        """
        # for finite volume we use the cell centres
        if symbol.name in ["x", "r"]:
            symbol_mesh = self.mesh.combine_submeshes(*symbol.domain)
            return pybamm.Vector(symbol_mesh[0].nodes)
        else:
            raise NotImplementedError("3D meshes not yet implemented")

    def broadcast(self, symbol, domain):
        """
        Broadcast symbol to a specified domain. To do this, calls
        :class:`pybamm.NumpyBroadcast`

        See :meth: `pybamm.SpatialMethod.broadcast`
        """
        return pybamm.NumpyBroadcast(symbol, domain, self.mesh)

    def gradient(self, symbol, discretised_symbol, boundary_conditions):
        """Matrix-vector multiplication to implement the gradient operator.
        See :meth:`pybamm.SpatialMethod.gradient`
        """
        # Check that boundary condition keys are hashes (ids)
        for key in boundary_conditions.keys():
            assert isinstance(key, int), TypeError(
                "boundary condition keys should be hashes, not {}".format(type(key))
            )
        # Discretise symbol
        domain = symbol.domain
        # Add Dirichlet boundary conditions, if defined
        if symbol.id in boundary_conditions:
            lbc = boundary_conditions[symbol.id]["left"]
            rbc = boundary_conditions[symbol.id]["right"]
            discretised_symbol = self.add_ghost_nodes(
                symbol, discretised_symbol, lbc, rbc
            )
            domain = (
                [domain[0] + "_left ghost cell"]
                + domain
                + [domain[-1] + "_right ghost cell"]
            )

        # note in 1D spherical grad and normal grad are the same
        gradient_matrix = self.gradient_matrix(domain)
        return gradient_matrix @ discretised_symbol

    def gradient_matrix(self, domain):
        """
        Gradient matrix for finite volumes in the appropriate domain.
        Equivalent to grad(y) = (y[1:] - y[:-1])/dx

        Parameters
        ----------
        domain : list
            The domain(s) in which to compute the gradient matrix

        Returns
        -------
        :class:`pybamm.Matrix`
            The (sparse) finite volume gradient matrix for the domain
        """
        # Create appropriate submesh by combining submeshes in domain
        submesh_list = self.mesh.combine_submeshes(*domain)

        # can just use 1st entry of list to obtain the point etc
        submesh = submesh_list[0]

        # Create 1D matrix using submesh
        n = submesh.npts
        e = 1 / submesh.d_nodes
        data = np.vstack(
            [np.concatenate([-e, np.array([0])]), np.concatenate([np.array([0]), e])]
        )
        diags = np.array([0, 1])
        sub_matrix = spdiags(data, diags, n - 1, n)

        # second dim length
        second_dim_len = len(submesh_list)

        # generate full matrix from the submatrix
        matrix = kron(eye(second_dim_len), sub_matrix)

        return pybamm.Matrix(matrix)

    def divergence(self, symbol, discretised_symbol, boundary_conditions):
        """Matrix-vector multiplication to implement the divergence operator.
        See :meth:`pybamm.SpatialMethod.divergence`
        """
        # Check that boundary condition keys are hashes (ids)
        for key in boundary_conditions.keys():
            assert isinstance(key, int), TypeError(
                "boundary condition keys should be hashes, not {}".format(type(key))
            )

        domain = symbol.domain
        submesh_list = self.mesh.combine_submeshes(*domain)
        # create a bc vector of length equal to the number variables
        # (only has non zero entries for neumann bcs)

        prim_dim = submesh_list[0].npts
        second_dim = len(submesh_list)
        total_pts = prim_dim * second_dim

        # Add Neumann boundary conditions if defined
        if symbol.id in boundary_conditions:
            # TODO:these are symbols so need to check them
            lbc = boundary_conditions[symbol.id]["left"]
            rbc = boundary_conditions[symbol.id]["right"]

            # doing via loop so that it is easier to implement x varing bcs
            bcs_symbol = pybamm.Vector(np.array([]))  # empty vector
            for i in range(len(submesh_list)):
                # only the interior equations:
                interior = pybamm.Vector(np.zeros(prim_dim - 2))
                left = -lbc / pybamm.Vector(np.array([submesh_list[i].d_edges[0]]))
                right = rbc / pybamm.Vector(np.array([submesh_list[i].d_edges[-1]]))
                bcs_symbol = pybamm.NumpyConcatenation(
                    bcs_symbol, left, interior, right
                )

            # now we must create a matrix of size (npts * (npts -1) )
            # this is a different size to the one created when we have
            # flux boundary conditions so need a flag
            divergence_matrix = self.divergence_matrix(domain, bc_type="neumann")

            # only need interior edges for spherical neumann
            edges = submesh_list[0].edges[1:-1]

        else:
            divergence_matrix = self.divergence_matrix(domain)
            bcs_vec = np.zeros(total_pts)
            bcs_symbol = pybamm.Vector(bcs_vec)
            # need all edges for spherical dirichlet
            edges = submesh_list[0].edges

        # check for particle domain
        if submesh_list[0].coord_sys == "spherical polar":

            # create np.array of repeated submesh[0].nodes
            r_numpy = np.kron(np.ones(second_dim), submesh_list[0].nodes)
            r_edges_numpy = np.kron(np.ones(second_dim), edges)

            r = pybamm.Vector(r_numpy)
            r_edges = pybamm.Vector(r_edges_numpy)

            # for clarity, we are implicitly multiplying the the lbc by r^2=0
            # and the rbc by r^2=1. But lbc is 0 so we don't need to do
            # any r_edges^2 operations on bcs_symbol
            out = (1 / (r ** 2)) * (
                divergence_matrix @ ((r_edges ** 2) * discretised_symbol) + bcs_symbol
            )
        else:
            out = divergence_matrix @ discretised_symbol + bcs_symbol

        return out

    def divergence_matrix(self, domain, bc_type="dirichlet"):
        """
        Divergence matrix for finite volumes in the appropriate domain.
        Equivalent to div(N) = (N[1:] - N[:-1])/dx

        Parameters
        ----------
        domain : list
            The domain(s) in which to compute the divergence matrix

        Returns
        -------
        :class:`pybamm.Matrix`
            The (sparse) finite volume divergence matrix for the domain
        """
        # Create appropriate submesh by combining submeshes in domain
        submesh_list = self.mesh.combine_submeshes(*domain)

        # can just use 1st entry of list to obtain the point etc
        submesh = submesh_list[0]
        e = 1 / submesh.d_edges

        # Create matrix using submesh
        n = submesh.npts + 1
        if bc_type == "dirichlet":
            data = np.vstack(
                [
                    np.concatenate([-e, np.array([0])]),
                    np.concatenate([np.array([0]), e]),
                ]
            )
            diags = np.array([0, 1])
            sub_matrix = spdiags(data, diags, n - 1, n)
        elif bc_type == "neumann":
            # we don't have to act on bc fluxes which are now in
            # the bc vector
            data = np.vstack([-e[1:], e[:-1]])
            diags = np.array([-1, 0])
            sub_matrix = spdiags(data, diags, n - 1, n - 2)
        else:
            raise NotImplementedError(
                "Can only process Neumann or Dirichlet boundary conditions"
            )

        # repeat matrix for each node in secondary dimensions
        second_dim_len = len(submesh_list)
        # generate full matrix from the submatrix
        matrix = kron(eye(second_dim_len), sub_matrix)
        return pybamm.Matrix(matrix)

    def integral(self, domain, symbol, discretised_symbol):
        """Vector-vector dot product to implement the integral operator.
        See :meth:`pybamm.BaseDiscretisation.integral`
        """
        # Calculate integration vector
        integration_vector = self.definite_integral_vector(domain)
        # Check for particle domain
        if ("negative particle" or "positive particle") in symbol.domain:
            submesh_list = self.mesh.combine_submeshes(*symbol.domain)
            second_dim = len(submesh_list)
            r_numpy = np.kron(np.ones(second_dim), submesh_list[0].nodes)
            r = pybamm.Vector(r_numpy)
            out = 2 * np.pi * integration_vector @ (discretised_symbol * r)
        else:
            out = integration_vector @ discretised_symbol
        out.domain = []
        return out

    def definite_integral_vector(self, domain):
        """
        Vector for finite-volume implementation of the definite integral

        .. math::
            I = \\int_{a}^{b}\\!f(s)\\,ds

        for where :math:`a` and :math:`b` are the left-hand and right-hand boundaries of
        the domain respectively

        Parameters
        ----------
        domain : list
            The domain(s) of integration

        Returns
        -------
        :class:`pybamm.Vector`
            The finite volume integral vector for the domain
        """
        # Create appropriate submesh by combining submeshes in domain
        submesh_list = self.mesh.combine_submeshes(*domain)

        # Create vector of ones using submesh
        vector = np.array([])
        for submesh in submesh_list:
            vector = np.append(vector, submesh.d_edges * np.ones_like(submesh.nodes))

        return pybamm.Vector(vector)

    def add_ghost_nodes(self, symbol, discretised_symbol, lbc, rbc):
        """
        Add Dirichlet boundary conditions via ghost nodes.

        For a boundary condition "y = a at the left-hand boundary",
        we concatenate a ghost node to the start of the vector y with value "2*a - y1"
        where y1 is the value of the first node.
        Similarly for the right-hand boundary condition.

        Currently, Dirichlet boundary conditions can only be applied on state
        variables (e.g. concentration, temperature), and not on expressions.
        To access the value of the first node (y1), we create a "first_node" object
        which is a StateVector whose y_slice is the start of the y_slice of
        discretised_symbol.
        Similarly, the last node is a StateVector whose y_slice is the end of the
        y_slice of discretised_symbol

        Parameters
        ----------
        discretised_symbol : :class:`pybamm.StateVector` (size n)
            The discretised variable (a state vector) to which to add ghost nodes
        lbc : :class:`pybamm.Scalar`
            Dirichlet bouncary condition on the left-hand side
        rbc : :class:`pybamm.Scalar`
            Dirichlet bouncary condition on the right-hand side

        Returns
        -------
        :class:`pybamm.Concatenation` (size n+2)
            Concatenation of the variable (a state vector) and ghost nodes

        """
        assert isinstance(discretised_symbol, pybamm.StateVector), NotImplementedError(
            """discretised_symbol must be a StateVector, not {}""".format(
                type(discretised_symbol)
            )
        )

        # determine the y_slice sizes
        y_slice_start = discretised_symbol.y_slice.start
        y_slice_stop = discretised_symbol.y_slice.stop

        y = np.arange(y_slice_start, y_slice_stop)

        # reshape y_slices into more helpful form
        submesh_list = self.mesh.combine_submeshes(*symbol.domain)
        if isinstance(submesh_list[0].npts, list):
            NotImplementedError("Can only take in 1D primary directions")
        # size = [submesh_list[0].npts, len(submesh_list)]
        size = [len(submesh_list), submesh_list[0].npts]
        y = np.reshape(y, size)
        y_left = y[:, 0]
        y_right = y[:, -1]

        new_discretised_symbol = pybamm.Vector(np.array([]))  # starts empty

        for i in range(len(submesh_list)):
            y_slice_start = y_left[i]
            y_slice_stop = y_right[i]

            # left ghost cell
            first_node = pybamm.StateVector(slice(y_slice_start, y_slice_start + 1))
            left_ghost_cell = 2 * lbc - first_node
            # middle symbol
            sub_disc_symbol = pybamm.StateVector(slice(y_slice_start, y_slice_stop + 1))
            # right ghost cell
            last_node = pybamm.StateVector(slice(y_slice_stop, y_slice_stop + 1))
            right_ghost_cell = 2 * rbc - last_node

            concatenated_sub_disc_symbol = pybamm.NumpyConcatenation(
                left_ghost_cell, sub_disc_symbol, right_ghost_cell
            )

            new_discretised_symbol = pybamm.NumpyConcatenation(
                new_discretised_symbol, concatenated_sub_disc_symbol
            )

        return new_discretised_symbol

    def boundary_value(self, discretised_symbol, side):
        """
        Uses linear extrapolation to get the boundary value of a variable in the
        Finite Volume Method.

        Parameters
        -----------
        discretised_symbol : :class:`pybamm.StateVector`
            The discretised variable from which to calculate the boundary value
        side : string
            Which side to take the boundary value on ("left" or "right")

        Returns
        -------
        :class:`pybamm.Symbol`
            The variable representing the boundary value.
        """

        def linear_extrapolation(array):
            """Linearly extrapolates an array"""
            if side == "left":
                return array[0] + (array[0] - array[1]) / 2
            elif side == "right":
                return array[-1] + (array[-1] - array[-2]) / 2

        return BoundaryValueEvaluated(discretised_symbol, linear_extrapolation)

    def mass_matrix(self, symbol, boundary_conditions):
        """
        Calculates the mass matrix for a spatial method.

        Parameters
        ----------
        symbol: :class:`pybamm.Variable`
            The variable corresponding to the equation for which we are
            calculating the mass matrix.
        boundary_conditions : dict
            The boundary conditions of the model
            ({symbol.id: {"left": left bc, "right": right bc}})

        Returns
        -------
        :class:`pybamm.Matrix`
            The (sparse) mass matrix for the spatial method.
        """
        # NOTE: for different spatial methods the matrix may need to be adjusted
        # to account for Dirichlet boundary conditions. Here, we just have that
        # the mass matrix is the identity.

        # Create appropriate submesh by combining submeshes in domain
        submesh = self.mesh.combine_submeshes(*symbol.domain)

        # Get number of points in primary dimension
        n = submesh[0].npts

        # Create mass matrix for primary dimension
        prim_mass = eye(n)

        # Get number of points in secondary dimension
        sec_pts = len(submesh)

        mass = kron(eye(sec_pts), prim_mass)
        return pybamm.Matrix(mass)

    #######################################################
    # Can probably be moved outside of the spatial method
    ######################################################

    def compute_diffusivity(self, discretised_symbol):
        """
        Compute the diffusivity at cell edges, based on the diffusivity at cell nodes.
        For now we just take the arithemtic mean, though it may be better to take the
        harmonic mean based on [1].

        [1] Recktenwald, Gerald. "The control-volume finite-difference approximation to
        the diffusion equation." (2012).

        Parameters
        ----------
        discretised_symbol : :class:`pybamm.Symbol`
            Symbol to be averaged. When evaluated, this symbol returns either a scalar
            or an array of shape (n,), where n is the number of points in the mesh for
            the symbol's domain (n = self.mesh[symbol.domain].npts)

        Returns
        -------
        :class:`pybamm.NodeToEdge`
            Averaged symbol. When evaluated, this returns either a scalar or an array of
            shape (n-1,) as appropriate.
        """

        def arithmetic_mean(array):
            """Calculate the arithemetic mean of an array"""
            return (array[1:] + array[:-1]) / 2

        return pybamm.NodeToEdge(discretised_symbol, arithmetic_mean)


class BoundaryValueEvaluated(pybamm.SpatialOperator):
    """A node in the expression tree representing a unary operator that evaluates the
    value of its child at a boundary.

    Parameters
    ----------
    child : :class:`Symbol`
        child node
    boundary_function : method
        the function used to calculate the boundary value

    **Extends:** :class:`pybamm.SpatialOperator`
    """

    def __init__(self, child, boundary_function):
        """ See :meth:`pybamm.UnaryOperator.__init__()`. """
        super().__init__(
            "boundary value ({})".format(boundary_function.__name__), child
        )
        self._boundary_function = boundary_function
        # Domain of BoundaryValue must be ([]) so that expressions can be formed
        # of boundary values of variables in different domains
        self.domain = []

    def evaluate(self, t=None, y=None):
        """ See :meth:`pybamm.Symbol.evaluate()`. """
        evaluated_child = self.children[0].evaluate(t, y)
        return self._boundary_function(evaluated_child)


class NodeToEdge(pybamm.SpatialOperator):
    """A node in the expression tree representing a unary operator that evaluates the
    value of its child at cell edges by averaging the value at cell nodes.

    Parameters
    ----------

    child : :class:`Symbol`
        child node
    node_to_edge_function : method
        the function used to average; only acts if the child evaluates to a
        one-dimensional numpy array

    **Extends:** :class:`pybamm.SpatialOperator`
    """

    def __init__(self, child, node_to_edge_function):
        """ See :meth:`pybamm.UnaryOperator.__init__()`. """
        super().__init__(
            "node to edge ({})".format(node_to_edge_function.__name__), child
        )
        self._node_to_edge_function = node_to_edge_function

    @property
    def node_to_edge_function(self):
        return self._node_to_edge_function

    def evaluate(self, t=None, y=None):
        """ See :meth:`pybamm.Symbol.evaluate()`. """
        evaluated_child = self.children[0].evaluate(t, y)
        # If the evaluated child is a numpy array of shape (n,), do the averaging
        # NOTE: Doing this check every time might be slow?
        if isinstance(evaluated_child, np.ndarray) and len(evaluated_child.shape) == 1:
            return self._node_to_edge_function(evaluated_child)
        # If not, no need to average
        else:
            return evaluated_child
