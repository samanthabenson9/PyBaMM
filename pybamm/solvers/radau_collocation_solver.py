#
# CasADi Solver class
#
import casadi
import pybamm
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq


class ImplicitRandauSolver(pybamm.BaseSolver):
    """Solve a discretised model, using CasADi.

    **Extends**: :class:`pybamm.BaseSolver`

    Parameters
    ----------
    method : str, optional
        The method to use for solving the system ('cvodes', for ODEs, or 'idas', for
        DAEs). Default is 'idas'.
    mode : str
            How to solve the model (default is "safe"):

            - "fast": perform direct integration, without accounting for events. \
            Recommended when simulating a drive cycle or other simulation where \
            no events should be triggered.
            - "safe": perform step-and-check integration in global steps of size \
            dt_max, checking whether events have been triggered. Recommended for \
            simulations of a full charge or discharge.
            - "old safe": perform step-and-check integration in steps of size dt \
            for each dt in t_eval, checking whether events have been triggered.
    rtol : float, optional
        The relative tolerance for the solver (default is 1e-6).
    atol : float, optional
        The absolute tolerance for the solver (default is 1e-6).
    root_method : str or pybamm algebraic solver class, optional
        The method to use to find initial conditions (for DAE solvers).
        If a solver class, must be an algebraic solver class.
        If "casadi",
        the solver uses casadi's Newton rootfinding algorithm to find initial
        conditions. Otherwise, the solver uses 'scipy.optimize.root' with method
        specified by 'root_method' (e.g. "lm", "hybr", ...)
    root_tol : float, optional
        The tolerance for root-finding. Default is 1e-6.
    max_step_decrease_counts : float, optional
        The maximum number of times step size can be decreased before an error is
        raised. Default is 5.
    dt_max : float, optional
        The maximum global step size (in seconds) used in "safe" mode. If None
        the default value corresponds to a non-dimensional time of 0.01
        (i.e. ``0.01 * model.timescale_eval``).
    extra_options_setup : dict, optional
        Any options to pass to the CasADi integrator when creating the integrator.
        Please consult `CasADi documentation <https://tinyurl.com/y5rk76os>`_ for
        details.
    extra_options_call : dict, optional
        Any options to pass to the CasADi integrator when calling the integrator.
        Please consult `CasADi documentation <https://tinyurl.com/y5rk76os>`_ for
        details.

    """

    def __init__(
        self,
        mode="safe",
        rtol=1e-6,
        atol=1e-6,
        root_method="casadi",
        root_tol=1e-6,
        max_step_decrease_count=5,
        dt_max=None,
        extra_options_setup=None,
        extra_options_call=None,
    ):
        super().__init__("problem dependent", rtol, atol, root_method, root_tol)
        if mode in ["safe", "fast", "old safe"]:
            self.mode = mode
        else:
            raise ValueError(
                """
                invalid mode '{}'. Must be either 'safe' or 'old safe', for solving
                with events, or 'fast', for solving quickly without events""".format(
                    mode
                )
            )
        self.max_step_decrease_count = max_step_decrease_count
        self.dt_max = dt_max

        self.extra_options_setup = extra_options_setup or {}
        self.extra_options_call = extra_options_call or {}

        self.name = "Implicit Radau Solver solver with '{}' mode".format(mode)

        # Need to update or remove this
        pybamm.citations.register("Andersson2019")

    def _integrate(self, model, t_eval, inputs=None):
        """
        Solve a DAE model defined by residuals with initial conditions y0.

        Parameters
        ----------
        model : :class:`pybamm.BaseModel`
            The model whose solution to calculate.
        t_eval : numeric type
            The times at which to compute the solution
        inputs : dict, optional
            Any external variables or input parameters to pass to the model when solving
        """
        inputs = inputs or {}
        # convert inputs to casadi format
        inputs = casadi.vertcat(*[x for x in inputs.values()])

        if self.mode == "fast":
            integrator = self.get_integrator(model, t_eval, inputs)
            solution = self._run_integrator(integrator, model, model.y0, inputs, t_eval)
            solution.termination = "final time"
            return solution
        elif not model.events:
            pybamm.logger.info("No events found, running fast mode")
            integrator = self.get_integrator(model, t_eval, inputs)
            solution = self._run_integrator(integrator, model, model.y0, inputs, t_eval)
            solution.termination = "final time"
            return solution
        elif self.mode == "safe":
            y0 = model.y0
            if isinstance(y0, casadi.DM):
                y0 = y0.full().flatten()
            # Step-and-check
            t = t_eval[0]
            t_f = t_eval[-1]
            init_event_signs = np.sign(
                np.concatenate(
                    [event(t, y0, inputs) for event in model.terminate_events_eval]
                )
            )
            pybamm.logger.info("Start solving {} with {}".format(model.name, self.name))

            # Initialize solution
            solution = pybamm.Solution(np.array([t]), y0[:, np.newaxis])
            solution.solve_time = 0

            # Try to integrate in global steps of size dt_max. Note: dt_max must
            # be at least as big as the the biggest step in t_eval (multiplied
            # by some tolerance, here 1.01) to avoid an empty integration window below
            if self.dt_max:
                # Non-dimensionalise provided dt_max
                dt_max = self.dt_max / model.timescale_eval
            else:
                dt_max = 0.01
            dt_eval_max = np.max(np.diff(t_eval)) * 1.01
            dt_max = np.max([dt_max, dt_eval_max])
            while t < t_f:
                # Step
                solved = False
                count = 0
                dt = dt_max
                while not solved:
                    # Get window of time to integrate over (so that we return
                    # all the points in t_eval, not just t and t+dt)
                    t_window = np.concatenate(
                        ([t], t_eval[(t_eval > t) & (t_eval < t + dt)])
                    )
                    # Sometimes near events the solver fails between two time
                    # points in t_eval (i.e. no points t < t_i < t+dt for t_i
                    # in t_eval), so we simply integrate from t to t+dt
                    if len(t_window) == 1:
                        t_window = np.array([t, t + dt])

                    integrator = self.get_integrator(model, t_window, inputs)
                    # Try to solve with the current global step, if it fails then
                    # halve the step size and try again.
                    try:
                        current_step_sol = self._run_integrator(
                            integrator, model, y0, inputs, t_window
                        )
                        solved = True
                    except pybamm.SolverError:
                        dt /= 2
                        # also reduce maximum step size for future global steps
                        dt_max = dt
                    count += 1
                    if count >= self.max_step_decrease_count:
                        raise pybamm.SolverError(
                            """
                            Maximum number of decreased steps occurred at t={}. Try
                            solving the model up to this time only or reducing dt_max.
                            """.format(
                                t
                            )
                        )
                # Check most recent y to see if any events have been crossed
                new_event_signs = np.sign(
                    np.concatenate(
                        [
                            event(t, current_step_sol.y[:, -1], inputs)
                            for event in model.terminate_events_eval
                        ]
                    )
                )
                # Exit loop if the sign of an event changes
                # Locate the event time using a root finding algorithm and
                # event state using interpolation. The solution is then truncated
                # so that only the times up to the event are returned
                if (new_event_signs != init_event_signs).any():
                    # get the index of the events that have been crossed
                    event_ind = np.where(new_event_signs != init_event_signs)[0]
                    active_events = [model.terminate_events_eval[i] for i in event_ind]

                    # create interpolant to evaluate y in the current integration
                    # window
                    y_sol = interp1d(current_step_sol.t, current_step_sol.y)

                    # loop over events to compute the time at which they were triggered
                    t_events = [None] * len(active_events)
                    for i, event in enumerate(active_events):

                        def event_fun(t):
                            return event(t, y_sol(t), inputs)

                        if np.isnan(event_fun(current_step_sol.t[-1])[0]):
                            # bracketed search fails if f(a) or f(b) is NaN, so we
                            # need to find the times for which we can evaluate the event
                            times = [
                                t
                                for t in current_step_sol.t
                                if event_fun(t)[0] == event_fun(t)[0]
                            ]
                        else:
                            times = current_step_sol.t
                        # skip if sign hasn't changed
                        if np.sign(event_fun(times[0])) != np.sign(
                            event_fun(times[-1])
                        ):
                            t_events[i] = brentq(
                                lambda t: event_fun(t), times[0], times[-1]
                            )
                        else:
                            t_events[i] = np.nan

                    # t_event is the earliest event triggered
                    t_event = np.nanmin(t_events)
                    y_event = y_sol(t_event)

                    # return truncated solution
                    t_truncated = current_step_sol.t[current_step_sol.t < t_event]
                    y_trunctaed = current_step_sol.y[:, 0 : len(t_truncated)]
                    truncated_step_sol = pybamm.Solution(t_truncated, y_trunctaed)
                    # assign temporary solve time
                    truncated_step_sol.solve_time = np.nan
                    # append solution from the current step to solution
                    solution.append(truncated_step_sol)

                    solution.termination = "event"
                    solution.t_event = t_event
                    solution.y_event = y_event
                    break
                else:
                    # assign temporary solve time
                    current_step_sol.solve_time = np.nan
                    # append solution from the current step to solution
                    solution.append(current_step_sol)
                    # update time
                    t = t_window[-1]
                    # update y0
                    y0 = solution.y[:, -1]
            return solution
        elif self.mode == "old safe":
            y0 = model.y0
            if isinstance(y0, casadi.DM):
                y0 = y0.full().flatten()
            # Step-and-check
            t = t_eval[0]
            init_event_signs = np.sign(
                np.concatenate(
                    [event(t, y0, inputs) for event in model.terminate_events_eval]
                )
            )
            pybamm.logger.info("Start solving {} with {}".format(model.name, self.name))

            # Initialize solution
            solution = pybamm.Solution(np.array([t]), y0[:, np.newaxis])
            solution.solve_time = 0
            for dt in np.diff(t_eval):
                # Step
                solved = False
                count = 0
                while not solved:
                    integrator = self.get_integrator(
                        model, np.array([t, t + dt]), inputs
                    )
                    # Try to solve with the current step, if it fails then halve the
                    # step size and try again. This will make solution.t slightly
                    # different to t_eval, but shouldn't matter too much as it should
                    # only happen near events.
                    try:
                        current_step_sol = self._run_integrator(
                            integrator, model, y0, inputs, np.array([t, t + dt])
                        )
                        solved = True
                    except pybamm.SolverError:
                        dt /= 2
                    count += 1
                    if count >= self.max_step_decrease_count:
                        raise pybamm.SolverError(
                            """
                            Maximum number of decreased steps occurred at t={}. Try
                            solving the model up to this time only.
                            """.format(
                                t
                            )
                        )
                # Check most recent y
                new_event_signs = np.sign(
                    np.concatenate(
                        [
                            event(t, current_step_sol.y[:, -1], inputs)
                            for event in model.terminate_events_eval
                        ]
                    )
                )
                # Exit loop if the sign of an event changes
                if (new_event_signs != init_event_signs).any():
                    solution.termination = "event"
                    solution.t_event = solution.t[-1]
                    solution.y_event = solution.y[:, -1]
                    break
                else:
                    # assign temporary solve time
                    current_step_sol.solve_time = np.nan
                    # append solution from the current step to solution
                    solution.append(current_step_sol)
                    # update time
                    t += dt
                    # update y0
                    y0 = solution.y[:, -1]
            return solution

    def get_integrator(self, model, t_eval, inputs):
        y0 = model.y0
        rhs = model.casadi_rhs
        algebraic = model.casadi_algebraic

        # set up
        pp = casadi.MX.sym("p", inputs.shape[0])

        # End time
        tf = t_eval[-1]
        # Number of finite elements
        n_eval = t_eval.shape[
            0
        ]  # siegeljb 6/3/2020, this is off by 1, need to fix eventually.
        # Size of the finite elements
        h_eval = tf / (n_eval - 1)  # this should work, but need to verify.

        
        h=min(1e-3,h_eval)
        n=int(t_eval[-1]//h+1)

        # Dimensions
        n_x = rhs(t_eval[0], y0, inputs).shape[0]
        n_p = inputs.shape[0]
        n_z = algebraic(t_eval[0], y0, inputs).shape[0]

        # Declare variables
        x = casadi.SX.sym("x", n_x)  # state
        p = casadi.SX.sym("u", n_p)  # control
        z = casadi.SX.sym("z", n_z)  # algeb state
        t = casadi.SX.sym("t")  # time.
# declare a time variable here... todo.
        yfull=casadi.vertcat(x, z)

        ode = rhs(
            t, yfull, p
        )  # vertcat(0.7*x[1]+sin(2.5*z[0]),1.4*x[0]+cos(2.5*z[0]))
        alg = algebraic(
            t, yfull, p
        )  # vertcat(z[0]**2+x[1]**2-1)

        fx = casadi.Function("fx", [x, z, p, t], [ode])
        fz = casadi.Function("fz", [x, z, p, t], [alg])

        

        # 0 = fA(x, z, p, t) at all times by means of the implicit function theorem,
        # implying in particular that ∂fA/∂z must be invertible.

        # Number of finite elements
        # n = 100

        # Size of the finite elements
        # h = tf / n

        # Degree of interpolating polynomial
        d = 5

        # Choose collocation points
        tau_root = [0] + casadi.collocation_points(d, "legendre")

        # Coefficients of the collocation equation
        C = np.zeros((d + 1, d + 1))

        # Coefficients of the continuity equation
        D = np.zeros(d + 1)

        # Dimensionless time inside one control interval
        tau = casadi.SX.sym("tau")

        # For all collocation points
        for j in range(d + 1):
            # Construct Lagrange polynomials to get the polynomial basis at the collocation point
            L = 1
            for r in range(d + 1):
                if r != j:
                    L *= (tau - tau_root[r]) / (tau_root[j] - tau_root[r])

            # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
            lfcn = casadi.Function("lfcn", [tau], [L])
            D[j] = lfcn(1.0)

            # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
            tfcn = casadi.Function("tfcn", [tau], [casadi.tangent(L, tau)])
            for r in range(d + 1):
                C[j, r] = tfcn(tau_root[r])

        # Total number of variables for one finite element
        X0 = casadi.MX.sym("X0", n_x)
        # Z0 = MX.sym('Z0',nz)
        Z0_guess = casadi.MX.sym("Z0", n_z)
        # XZ0=MX.sym('XZ0',nx+nz)
        # X0=XZ0[0:nx]
        # Z0=XZ0[nx:(nx+nz)]
        T = casadi.MX.sym("T")
        P = casadi.MX.sym("P", n_p)
        V = casadi.MX.sym("V", (d * n_x + (d + 1) * n_z))
        Vx = V[0 : d * n_x]  # MX.sym('Vz',d*nx)

        if algebraic(t_eval[0], y0, pp).is_empty():

            



            # # End time
            # tf = 1.0

            # # Dimensions
            # n_x = 1
            # n_p = 0

            # # Declare variables
            # x = casadi.SX.sym("x", n_x)  # state
            # p = casadi.SX.sym("u", n_p)  # control

            # ODE right hand side function
            # ode = vertcat((1 - x[1]*x[1])*x[0] - x[1] + p, \
            #   x[0], \
            #   x[0]*x[0] + x[1]*x[1] + p*p)

            # ode = 0.1 * x[0]
            # dae = {"x": x, "p": p, "ode": ode}
            # f = casadi.Function("f", [x, p], [ode])

            f = casadi.Function("f", [x, p, t], [ode])
         

            # # Number of finite elements
            # n = 100

            # # Size of the finite elements
            # h = tf / n

            # # Degree of interpolating polynomial
            # d = 4

            # # Choose collocation points
            # tau_root = [0] + casadi.collocation_points(d, "legendre")

            # # Coefficients of the collocation equation
            # C = np.zeros((d + 1, d + 1))

            # # Coefficients of the continuity equation
            # D = np.zeros(d + 1)

            # # Dimensionless time inside one control interval
            # tau = casadi.SX.sym("tau")

            # # For all collocation points
            # for j in range(d + 1):
            #     # Construct Lagrange polynomials to get the polynomial basis at the collocation point
            #     L = 1
            #     for r in range(d + 1):
            #         if r != j:
            #             L *= (tau - tau_root[r]) / (tau_root[j] - tau_root[r])

            #     # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
            #     lfcn = casadi.Function("lfcn", [tau], [L])
            #     D[j] = lfcn(1.0)

            #     # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
            #     tfcn = casadi.Function("tfcn", [tau], [casadi.tangent(L, tau)])
            #     for r in range(d + 1):
            #         C[j, r] = tfcn(tau_root[r])

            # Total number of variables for one finite element
            X0 = casadi.MX.sym("X0", n_x)
            P = casadi.MX.sym("P", n_p)
            T = casadi.MX.sym("T")
            V = casadi.MX.sym("V", d * n_x)

            # Get the state at each collocation point
            X = [X0] + casadi.vertsplit(V, [r * n_x for r in range(d + 1)])

            # Get the collocation quations (that define V)
            V_eq = []
            for j in range(1, d + 1):
                # Expression for the state derivative at the collocation point
                xp_j = 0
                for r in range(d + 1):
                    xp_j += C[r, j] * X[r]

                # Append collocation equations
                f_j = f(X[j], P,T)
                V_eq.append(h * f_j - xp_j)

            # Concatenate constraints
            V_eq = casadi.vertcat(*V_eq)

            # Root-finding function, implicitly defines V as a function of X0 and P
            vfcn = casadi.Function("vfcn", [V, X0, P,T], [V_eq])

            # Convert to SX to decrease overhead
            vfcn_sx = vfcn.expand()

            # Create a implicit function instance to solve the system of equations
            # ifcn = casadi.rootfinder("ifcn", "fast_newton", vfcn_sx)
            opts = {}
            opts["max_iter"] = 5000
            # opts['reltol']=1e-6
            opts["abstol"] = 1e-8
            # Create a implicit function instance to solve the system of equations
            ifcn = casadi.rootfinder("ifcn", "fast_newton", vfcn_sx, opts)

            V = ifcn(casadi.MX(), X0, P,T)
            X = [X0 if r == 0 else V[(r - 1) * n_x : r * n_x] for r in range(d + 1)]

            # Get an expression for the state at the end of the finie element
            XF = 0
            for r in range(d + 1):
                XF += D[r] * X[r]

            # Get the discrete time dynamics
            F = casadi.Function("F", [X0, P, T], [XF])

            # Do this iteratively for all finite elements
            # Xs = casadi.MX(n_x, n)
            # Zs = casadi.MX(n_z, n)
            Xs = casadi.MX(n_x, n_eval)
            Zs = casadi.MX(n_z, n_eval)
            XsP = casadi.MX(n_x, 1)
            ZsP = casadi.MX(n_z, 1)
                

            Xs[:, 0] = X0
            X = X0
            # ts== casadi.MX(1 n)
            # ts[0]=0
            k=0
            for i in range(0,n-1):
                XsP=X
                # X = F(X, P,t_eval[i])
                X = F(X, P,i*h) #t_eval[k]
                if(i*h>=k*h_eval):
                    k=k+1
                    Xs[:, k] = X # should take the weighted average here, rather than the next point TODO siegeljb...
                #ts[i+1]
            # Take last step and update output
            X = F(X, P,tf)
            Xs[:, -1] = X

            # Fixed-step integrator
            irk_integrator = casadi.Function(
                "irk_integrator", {"x0": X0,'z0':Z0_guess,"p": P, "xf": Xs, "zf": Zs}, casadi.integrator_in(), casadi.integrator_out()
            )
            # sol_irk =irk_integrator(x0=1, p=[])
            # print(sol_irk['xf'])
            # irk_integrator = casadi.Function(
            #     "irk_integrator",
            #     {"x0": X0, "z0": Z0_guess, "p": P, "xf": Xs, "zf": Zs},
            #     casadi.integrator_in(),
            #     casadi.integrator_out(),
            # )

        else:

            # add initial guess for z0 here..
            Vz = V[d * n_x : (d * n_x + (d + 1) * n_z)]  # MX.sym('Vz',nz)

            # Get the state at each collocation point
            X = [X0] + casadi.vertsplit(Vx, [r * n_x for r in range(d + 1)])
            # Z= [Z0] + vertsplit(Vz,[r*nz for r in range(d+1)])
            Z = casadi.vertsplit(Vz, [r * n_z for r in range(d + 2)])

            #
            # if nz == 1:
            #   Z = [Z0] + [Vz]
            # else:
            #   Z= [Z0] + Vz

            # Get the collocation equations (that define V)
            V_eq = []
            for j in range(1, d + 1):
                # Expression for the state derivative at the collocation point
                xp_j = 0
                #  zp_j = 0
                for r in range(d + 1):
                    xp_j += C[r, j] * X[r]
                #  zp_j += C[r,j]*Z[r]
                # Append collocation equations
                f_j = fx(X[j], Z[j], P,T)
                V_eq.append(h * f_j - xp_j)

            for j in range(d + 1):
                # Append algebraic constraints
                # fz_j=fz(X[j],Z[j],P)
                V_eq.append(fz(X[j], Z[j], P,T))

            # add inital constraint for Z0
            # V_eq.append(fz(X0, Z0, P))
            # Concatenate constraints
            V_eq = casadi.vertcat(*V_eq)

            # Root-finding function, implicitly defines V as a function of X0 and P

            vfcn = casadi.Function("vfcn", [V, X0, P, T], [V_eq])
            # vfcn = Function('vfcn', [V, X0,Z0_guess,P], [V_eq])

            # Convert to SX to decrease overhead
            vfcn_sx = vfcn.expand()

            opts = {}
            opts["max_iter"] = 5000
            # opts['reltol']=1e-6
            opts["abstol"] = 1e-8
            # Create a implicit function instance to solve the system of equations
            ifcn = casadi.rootfinder("ifcn", "fast_newton", vfcn_sx, opts)
            init_guess = casadi.MX.zeros(d * n_x + (d + 1) * n_z, 1)
            init_guess[d * n_x : (d * n_x + n_z)] = Z0_guess
            V = ifcn(init_guess, X0, P,T)
            # X = [X0 if r==0 else V[(r-1)*nx:r*nx] for r in range(d+1)]

            X = [X0 if r == 0 else V[(r - 1) * (n_x) : r * n_x] for r in range(d + 1)]
            # Z= [Z0 if r==0 else V[(d*nx+(r-1)*(nz)):(d*nx+r*nz)] for r in range(d+1)]
            Z = [
                V[(d * n_x + (r) * (n_z)) : (d * n_x + (r + 1) * n_z)]
                for r in range(d + 1)
            ]
            Z0 = Z[0]
            # Get an expression for the state at the end of the finie element
            XF = 0
            ZF = 0
            for r in range(d + 1):
                XF += D[r] * X[r]
                ZF += D[r] * Z[r]

            # Get the discrete time dynamics
            # F = Function('F', [X0,P],[XF,ZF,Z0])
            F = casadi.Function(
                "F", [X0, Z0_guess, P,T], [XF, ZF, Z0] #
            )  #

            # Do this iteratively for all finite elements
            # # change back for symbolic eval of integrator
            Xs = casadi.MX(n_x, n_eval)
            Zs = casadi.MX(n_z, n_eval)
            Xs[:, 0] = X0
            Zs[:, 0] = Z0_guess
            Xss=X0
            Zss= Z0_guess
            XZ=F(Xss, Zss, P,0)
            Zs[:, 0] = XZ[2]  # update inital guess for algebraic state.
            Zss=Zs[:, 0]
            k=0
            for i in range(n - 1):
                XZprev=XZ
                XZ = F(Xss, Zss, P,i*h)    

                Xss = XZ[0]
                Zss = XZ[1]
                if(i*h>=k*h_eval):
                    k=k+1
                    Xs[:, k] = Xss
                    Zs[:, k] = Zss
            # and return all values
            XZ = F(Xss, Zss, P,n*h)
            Xss = XZ[0]
            Zss = XZ[1]
            Xs[:, -1] = Xss
            Zs[:, -1] = Zss
            # # Test values
            # x0_val = np.array([1, 0])

            # # x0_val  = np.array([1,0])
            # # x0_val  = np.array([3])
            # z0_val = np.array([0])

            # p_val = 0.2

            # ts = np.linspace(0, tf, n)
            # #integrator = integrator('integrator', 'cvodes', dae, {'grid':ts, 'output_t0':True})
            # mirk_integrator = Function('mirk_integrator', {'x0':X0, 'z0':Z0, 'p':P, 'xf':X},
            #                           integrator_in(), integrator_out())

            irk_integrator = casadi.Function(
                "irk_integrator",
                {"x0": X0, "z0": Z0_guess, "p": P, "xf": Xs, "zf": Zs},
                casadi.integrator_in(),
                casadi.integrator_out(),
            )

        return irk_integrator

    def _run_integrator(self, integrator, model, y0, inputs, t_eval):
        rhs_size = model.concatenated_rhs.size
        y0_diff, y0_alg = np.split(y0, [rhs_size])
        try:
            # Try solving
            sol = integrator(x0=y0_diff, z0=y0_alg, p=inputs, **self.extra_options_call)

            y_values = np.concatenate([sol["xf"].full(), sol["zf"].full()])
            return pybamm.Solution(t_eval, y_values)

        except RuntimeError as e:
            # If it doesn't work raise error
            raise pybamm.SolverError(e.args[0])


#
# sol_irk = irk_integrator(x0=x_sim, p=[0], z0=z_sim)


# plt.subplot(211)
# plt.plot(ts, N.array(sol['xf'])[0,:],label='idas')
# plt.plot(ts, N.array(sol_irk['xf'])[0,:],label='colloc')
# plt.legend()

# plt.subplot(212)
# plt.plot(ts, N.array(sol['zf'])[0,:],label='idas')
# plt.plot(ts, N.array(sol_irk['zf'])[0,:],label='colloc')
# plt.legend()
# plt.xlabel('t')
# plt.ylabel('z')
# plt.show()

# a=0

#! Plot the solution

# plt.plot(ts, N.array(sol["xf"])[0, :], label="idas")
# plt.plot(ts, N.array(sol["xf"])[1, :], label="idas")
# plt.plot(ts, N.array(sol_irk["xf"])[0, :], label="colloc")
# plt.plot(ts, N.array(sol_irk["xf"])[1, :], label="colloc")
# plt.plot(ts, N.array(sol["zf"])[0, :], label="idas")
# plt.plot(ts, N.array(sol_irk["zf"])[0, :], label="colloc")
# plt.legend()
# plt.xlabel("t")
# plt.ylabel("x")
# plt.show()
