#
# Simulations: discharge of a lead-acid battery
#
import pybamm


def options_to_tuple(options):
    bc_options = tuple(options["bc_options"].items())
    other_options = tuple(
        {k: v for k, v in options.items() if k != "bc_options"}.items()
    )
    return (*bc_options, *other_options)


def model_comparison(models, Crates, t_eval):
    " Solve models at a range of Crates "
    # load parameter values and geometry
    geometry = models[0].default_geometry
    param = models[0].default_parameter_values

    # Process parameters (same parameters for all models)
    for model in models:
        param.process_model(model)
    param.process_geometry(geometry)

    # set mesh
    var = pybamm.standard_spatial_vars
    var_pts = {var.x_n: 200, var.x_s: 200, var.x_p: 200}
    mesh = pybamm.Mesh(geometry, models[-1].default_submesh_types, var_pts)

    # discretise models
    discs = {}
    for model in models:
        disc = pybamm.Discretisation(mesh, model.default_spatial_methods)
        disc.process_model(model)
        # Store discretisation
        discs[model] = disc

    # solve model for range of Crates
    all_variables = {}
    for Crate in Crates:
        all_variables[Crate] = {}
        current = Crate * 17
        pybamm.logger.info("Setting typical current to {} A".format(current))
        param.update({"Typical current [A]": current})
        for model in models:
            param.update_model(model, discs[model])
            solution = model.default_solver.solve(model, t_eval)
            vars = pybamm.post_process_variables(
                model.variables, solution.t, solution.y, mesh
            )
            vars["solution"] = solution
            all_variables[Crate][(model.name, options_to_tuple(model.options))] = vars

    return all_variables, t_eval


def convergence_study(models, Crate, t_eval, all_npts):
    " Solve models at a range of number of grid points "
    # load parameter values and geometry
    geometry = models[0].default_geometry
    param = models[0].default_parameter_values
    current = Crate * 17
    param.update({"Typical current [A]": current})

    # Process parameters (same parameters for all models)
    for model in models:
        param.process_model(model)
    param.process_geometry(geometry)

    # set mesh
    var = pybamm.standard_spatial_vars

    # solve model for range of Crates
    all_variables = {}
    for npts in all_npts:
        all_variables[npts] = {}
        pybamm.logger.info("Setting number of grid points to {}".format(npts))
        var_pts = {var.x_n: npts, var.x_s: npts, var.x_p: npts}
        mesh = pybamm.Mesh(geometry, models[-1].default_submesh_types, var_pts)
        disc = pybamm.Discretisation(mesh, model.default_spatial_methods)

        # discretise models
        discs = {}
        for model in models:
            model_disc = disc.process_model(model, inplace=False)
            solution = model.default_solver.solve(model_disc, t_eval)
            vars = pybamm.post_process_variables(
                model_disc.variables, solution.t, solution.y, mesh
            )
            vars["solution"] = solution
            all_variables[npts][(model.name, options_to_tuple(model.options))] = vars

    return all_variables, t_eval
