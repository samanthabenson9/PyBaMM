#
# Test the base experiment class
#
import pybamm
import numpy as np
import unittest
import pandas as pd
import os


class TestExperiment(unittest.TestCase):
    def test_read_strings(self):
        # Import drive cycle from file
        drive_cycle = pd.read_csv(
            pybamm.get_parameters_filepath(
                os.path.join("input", "drive_cycles", "US06.csv")
            ),
            comment="#",
            header=None,
        ).to_numpy()

        experiment = pybamm.Experiment(
            [
                "Discharge at 1C for 0.5 hours",
                "Discharge at C/20 for 0.5 hours",
                "Charge at 0.5 C for 45 minutes",
                "Discharge at 1 A for 0.5 hours",
                "Charge at 200 mA for 45 minutes (1 minute period)",
                "Discharge at 1W for 0.5 hours",
                "Charge at 200mW for 45 minutes",
                "Rest for 10 minutes (5 minute period)",
                "Hold at 1V for 20 seconds",
                "Charge at 1 C until 4.1V",
                "Hold at 4.1 V until 50mA",
                "Hold at 3V until C/50",
                "Discharge at C/3 for 2 hours or until 2.5 V",
                "Run US06 (A)",
                "Run US06 (V) for 5 minutes",
                "Run US06 (W) for 0.5 hours",
            ],
            drive_cycles={"US06": drive_cycle},
            period="20 seconds",
        )

        self.assertEqual(
            experiment.operating_conditions[:-3],
            [
                {
                    "C-rate input [-]": 1,
                    "type": "C-rate",
                    "time": 1800.0,
                    "period": 20.0,
                    "dc_data": None,
                    "string": "Discharge at 1C for 0.5 hours",
                    "events": None,
                },
                {
                    "C-rate input [-]": 0.05,
                    "type": "C-rate",
                    "time": 1800.0,
                    "period": 20.0,
                    "dc_data": None,
                    "string": "Discharge at C/20 for 0.5 hours",
                    "events": None,
                },
                {
                    "C-rate input [-]": -0.5,
                    "type": "C-rate",
                    "time": 2700.0,
                    "period": 20.0,
                    "dc_data": None,
                    "string": "Charge at 0.5 C for 45 minutes",
                    "events": None,
                },
                {
                    "Current input [A]": 1,
                    "type": "current",
                    "time": 1800.0,
                    "period": 20.0,
                    "dc_data": None,
                    "string": "Discharge at 1 A for 0.5 hours",
                    "events": None,
                },
                {
                    "Current input [A]": -0.2,
                    "type": "current",
                    "time": 2700.0,
                    "period": 60.0,
                    "dc_data": None,
                    "string": "Charge at 200 mA for 45 minutes",
                    "events": None,
                },
                {
                    "Power input [W]": 1,
                    "type": "power",
                    "time": 1800.0,
                    "period": 20.0,
                    "dc_data": None,
                    "string": "Discharge at 1W for 0.5 hours",
                    "events": None,
                },
                {
                    "Power input [W]": -0.2,
                    "type": "power",
                    "time": 2700.0,
                    "period": 20.0,
                    "dc_data": None,
                    "string": "Charge at 200mW for 45 minutes",
                    "events": None,
                },
                {
                    "Current input [A]": 0,
                    "type": "current",
                    "time": 600.0,
                    "period": 300.0,
                    "dc_data": None,
                    "string": "Rest for 10 minutes",
                    "events": None,
                },
                {
                    "Voltage input [V]": 1,
                    "type": "voltage",
                    "time": 20.0,
                    "period": 20.0,
                    "dc_data": None,
                    "string": "Hold at 1V for 20 seconds",
                    "events": None,
                },
                {
                    "C-rate input [-]": -1,
                    "type": "C-rate",
                    "time": None,
                    "period": 20.0,
                    "dc_data": None,
                    "string": "Charge at 1 C until 4.1V",
                    "events": {"Voltage input [V]": 4.1, "type": "voltage"},
                },
                {
                    "Voltage input [V]": 4.1,
                    "type": "voltage",
                    "time": None,
                    "period": 20.0,
                    "dc_data": None,
                    "string": "Hold at 4.1 V until 50mA",
                    "events": {"Current input [A]": 0.05, "type": "current"},
                },
                {
                    "Voltage input [V]": 3,
                    "type": "voltage",
                    "time": None,
                    "period": 20.0,
                    "dc_data": None,
                    "string": "Hold at 3V until C/50",
                    "events": {"C-rate input [-]": 0.02, "type": "C-rate"},
                },
                {
                    "C-rate input [-]": 1 / 3,
                    "type": "C-rate",
                    "time": 7200.0,
                    "period": 20.0,
                    "dc_data": None,
                    "string": "Discharge at C/3 for 2 hours or until 2.5 V",
                    "events": {"Voltage input [V]": 2.5, "type": "voltage"},
                },
            ],
        )
        # Calculation for operating conditions of drive cycle
        time_0 = drive_cycle[:, 0][-1]
        period_0 = np.min(np.diff(drive_cycle[:, 0]))
        drive_cycle_1 = experiment.extend_drive_cycle(drive_cycle, end_time=300)
        time_1 = drive_cycle_1[:, 0][-1]
        period_1 = np.min(np.diff(drive_cycle_1[:, 0]))
        drive_cycle_2 = experiment.extend_drive_cycle(drive_cycle, end_time=1800)
        time_2 = drive_cycle_2[:, 0][-1]
        period_2 = np.min(np.diff(drive_cycle_2[:, 0]))
        # Check drive cycle operating conditions
        np.testing.assert_array_equal(
            experiment.operating_conditions[-3]["dc_data"], drive_cycle
        )
        self.assertEqual(experiment.operating_conditions[-3]["type"], "current")
        self.assertEqual(experiment.operating_conditions[-3]["time"], time_0)
        self.assertEqual(experiment.operating_conditions[-3]["period"], period_0)
        np.testing.assert_array_equal(
            experiment.operating_conditions[-2]["dc_data"], drive_cycle_1
        )
        self.assertEqual(experiment.operating_conditions[-2]["type"], "voltage")
        self.assertEqual(experiment.operating_conditions[-2]["time"], time_1)
        self.assertEqual(experiment.operating_conditions[-2]["period"], period_1)
        np.testing.assert_array_equal(
            experiment.operating_conditions[-1]["dc_data"], drive_cycle_2
        )
        self.assertEqual(experiment.operating_conditions[-1]["type"], "power")
        self.assertEqual(experiment.operating_conditions[-1]["time"], time_2)
        self.assertEqual(experiment.operating_conditions[-1]["period"], period_2)
        self.assertEqual(experiment.period, 20)

    def test_read_strings_cccv_combined(self):
        experiment = pybamm.Experiment(
            [
                (
                    "Discharge at C/20 for 0.5 hours",
                    "Charge at 0.5 C until 1V",
                    "Hold at 1V until C/50",
                    "Discharge at C/20 for 0.5 hours",
                ),
            ],
            cccv_handling="ode",
        )
        self.assertEqual(
            experiment.operating_conditions,
            [
                {
                    "C-rate input [-]": 0.05,
                    "type": "C-rate",
                    "time": 1800.0,
                    "period": 60.0,
                    "dc_data": None,
                    "string": "Discharge at C/20 for 0.5 hours",
                    "events": None,
                },
                {
                    "type": "CCCV",
                    "C-rate input [-]": -0.5,
                    "Voltage input [V]": 1,
                    "time": None,
                    "period": 60.0,
                    "dc_data": None,
                    "string": "Charge at 0.5 C until 1V then hold at 1V until C/50",
                    "events": {"C-rate input [-]": 0.02, "type": "C-rate"},
                },
                {
                    "C-rate input [-]": 0.05,
                    "type": "C-rate",
                    "time": 1800.0,
                    "period": 60.0,
                    "dc_data": None,
                    "string": "Discharge at C/20 for 0.5 hours",
                    "events": None,
                },
            ],
        )

        # Cases that don't quite match shouldn't do CCCV setup
        experiment = pybamm.Experiment(
            [
                "Charge at 0.5 C until 2V",
                "Hold at 1V until C/50",
            ],
            cccv_handling="ode",
        )
        self.assertEqual(
            experiment.operating_conditions,
            [
                {
                    "C-rate input [-]": -0.5,
                    "type": "C-rate",
                    "time": None,
                    "period": 60.0,
                    "dc_data": None,
                    "string": "Charge at 0.5 C until 2V",
                    "events": {"Voltage input [V]": 2, "type": "voltage"},
                },
                {
                    "Voltage input [V]": 1,
                    "type": "voltage",
                    "time": None,
                    "period": 60.0,
                    "dc_data": None,
                    "string": "Hold at 1V until C/50",
                    "events": {"C-rate input [-]": 0.02, "type": "C-rate"},
                },
            ],
        )
        experiment = pybamm.Experiment(
            [
                "Charge at 0.5 C for 2 minutes",
                "Hold at 1V until C/50",
            ],
            cccv_handling="ode",
        )
        self.assertEqual(
            experiment.operating_conditions,
            [
                {
                    "C-rate input [-]": -0.5,
                    "type": "C-rate",
                    "time": 120.0,
                    "period": 60.0,
                    "dc_data": None,
                    "string": "Charge at 0.5 C for 2 minutes",
                    "events": None,
                },
                {
                    "Voltage input [V]": 1,
                    "type": "voltage",
                    "time": None,
                    "period": 60.0,
                    "dc_data": None,
                    "string": "Hold at 1V until C/50",
                    "events": {"C-rate input [-]": 0.02, "type": "C-rate"},
                },
            ],
        )

    def test_cycle_unpacking(self):
        experiment = pybamm.Experiment(
            [
                ("Discharge at C/20 for 0.5 hours", "Charge at C/5 for 45 minutes"),
                ("Discharge at C/20 for 0.5 hours"),
                "Charge at C/5 for 45 minutes",
            ]
        )
        self.assertEqual(
            experiment.operating_conditions,
            [
                {
                    "C-rate input [-]": 0.05,
                    "type": "C-rate",
                    "time": 1800.0,
                    "period": 60.0,
                    "dc_data": None,
                    "string": "Discharge at C/20 for 0.5 hours",
                    "events": None,
                },
                {
                    "C-rate input [-]": -0.2,
                    "type": "C-rate",
                    "time": 2700.0,
                    "period": 60.0,
                    "dc_data": None,
                    "string": "Charge at C/5 for 45 minutes",
                    "events": None,
                },
                {
                    "C-rate input [-]": 0.05,
                    "type": "C-rate",
                    "time": 1800.0,
                    "period": 60.0,
                    "dc_data": None,
                    "string": "Discharge at C/20 for 0.5 hours",
                    "events": None,
                },
                {
                    "C-rate input [-]": -0.2,
                    "type": "C-rate",
                    "time": 2700.0,
                    "period": 60.0,
                    "dc_data": None,
                    "string": "Charge at C/5 for 45 minutes",
                    "events": None,
                },
            ],
        )
        self.assertEqual(experiment.cycle_lengths, [2, 1, 1])

    def test_str_repr(self):
        conds = ["Discharge at 1 C for 20 seconds", "Charge at 0.5 W for 10 minutes"]
        experiment = pybamm.Experiment(conds)
        self.assertEqual(
            str(experiment),
            "[('Discharge at 1 C for 20 seconds',)"
            + ", ('Charge at 0.5 W for 10 minutes',)]",
        )
        self.assertEqual(
            repr(experiment),
            "pybamm.Experiment([('Discharge at 1 C for 20 seconds',)"
            + ", ('Charge at 0.5 W for 10 minutes',)])",
        )

    def test_bad_strings(self):
        with self.assertRaisesRegex(ValueError, "cccv_handling"):
            pybamm.Experiment(
                ["Discharge at 1 C for 20 seconds", "Charge at 0.5 W for 10 minutes"],
                cccv_handling="bad",
            )
        with self.assertRaisesRegex(
            TypeError, "Operating conditions should be strings or tuples of strings"
        ):
            pybamm.Experiment([1, 2, 3])
        with self.assertRaisesRegex(
            TypeError, "Operating conditions should be strings or tuples of strings"
        ):
            pybamm.Experiment([(1, 2, 3)])
        with self.assertRaisesRegex(ValueError, "Operating conditions must contain"):
            pybamm.Experiment(["Discharge at 1 A at 2 hours"])
        with self.assertRaisesRegex(ValueError, "Instruction must be"):
            pybamm.Experiment(["Run at 1 A for 2 hours"])
        with self.assertRaisesRegex(ValueError, "Type of drive cycle must be"):
            pybamm.Experiment(["Run US06 for 2 hours"])
        with self.assertRaisesRegex(ValueError, "Instruction must be"):
            pybamm.Experiment(["Run at at 1 A for 2 hours"])
        with self.assertRaisesRegex(ValueError, "Instruction must be"):
            pybamm.Experiment(["Play at 1 A for 2 hours"])
        with self.assertRaisesRegex(ValueError, "Instruction"):
            pybamm.Experiment(["Cell Charge at 1 A for 2 hours"])
        with self.assertRaisesRegex(ValueError, "units must be"):
            pybamm.Experiment(["Discharge at 1 B for 2 hours"])
        with self.assertRaisesRegex(ValueError, "time units must be"):
            pybamm.Experiment(["Discharge at 1 A for 2 years"])

    def test_termination(self):
        experiment = pybamm.Experiment(["Discharge at 1 C for 20 seconds"])
        self.assertEqual(experiment.termination, {})

        experiment = pybamm.Experiment(
            ["Discharge at 1 C for 20 seconds"], termination="80.7% capacity"
        )
        self.assertEqual(experiment.termination, {"capacity": (80.7, "%")})

        experiment = pybamm.Experiment(
            ["Discharge at 1 C for 20 seconds"], termination="80.7 % capacity"
        )
        self.assertEqual(experiment.termination, {"capacity": (80.7, "%")})

        experiment = pybamm.Experiment(
            ["Discharge at 1 C for 20 seconds"], termination="4.1Ah capacity"
        )
        self.assertEqual(experiment.termination, {"capacity": (4.1, "Ah")})

        experiment = pybamm.Experiment(
            ["Discharge at 1 C for 20 seconds"], termination="4.1 A.h capacity"
        )
        self.assertEqual(experiment.termination, {"capacity": (4.1, "Ah")})

        experiment = pybamm.Experiment(
            ["Discharge at 1 C for 20 seconds"], termination="3V"
        )
        self.assertEqual(experiment.termination, {"voltage": (3, "V")})

        experiment = pybamm.Experiment(
            ["Discharge at 1 C for 20 seconds"], termination=["3V", "4.1Ah capacity"]
        )
        self.assertEqual(
            experiment.termination, {"voltage": (3, "V"), "capacity": (4.1, "Ah")}
        )

        with self.assertRaisesRegex(ValueError, "Only capacity"):
            experiment = pybamm.Experiment(
                ["Discharge at 1 C for 20 seconds"], termination="bla bla capacity bla"
            )
        with self.assertRaisesRegex(ValueError, "Only capacity"):
            experiment = pybamm.Experiment(
                ["Discharge at 1 C for 20 seconds"], termination="4 A.h something else"
            )
        with self.assertRaisesRegex(ValueError, "Capacity termination"):
            experiment = pybamm.Experiment(
                ["Discharge at 1 C for 20 seconds"], termination="1 capacity"
            )


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    pybamm.settings.debug_mode = True
    unittest.main()
