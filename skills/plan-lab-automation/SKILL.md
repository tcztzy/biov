---
name: plan-lab-automation
description: Review and simulate PyLabRobot scripts, labware layouts, tip use, and liquid tracking. Use for preparing or checking laboratory automation code without operating physical equipment.
---

# Review and simulate laboratory automation

Use the selected PyLabRobot version's
[official documentation](https://docs.pylabrobot.org/stable/) for liquid handling,
resource definitions, deck coordinates and material movement. The upstream
bundled tutorial was a snapshot; consult the current API for the pinned version.

1. Read the complete script and its imported local modules before execution.
   Identify every backend and side effect. Syntax-checking an AST does not
   sandbox code, and a timeout thread does not stop a device command.
2. Prepare a separate simulation script that explicitly constructs
   `LiquidHandlerChatterboxBackend` (or a documented simulation backend for that
   device). Use the actual supplied labware/deck definitions. Do not patch the
   text `STARBackend()` and assume other hardware connections are disabled.
3. Enable native tip and volume tracking; initialize known starting tips,
   liquid volumes, capacities and deck positions. Unknown inventory is missing
   input, not an empty/default reservoir. Run only the reviewed simulation in a
   separate process with a bounded runtime and checked exit status.
4. Inspect transfers for tip availability, source/destination capacity,
   aspiration/dispense ordering, compatible channels and resource collisions
   covered by the simulator. Save the simulator's actual commands and tracking
   results; do not invent operation/volume counts from a success message.

Return a JSON summary with `mode: "simulation"`, script and dependency versions,
configured backend, initial inventory, checks performed, observed failures and
output paths. Simulation verifies software/resource assumptions, not physical
calibration, sterility or liquid-class performance. Physical execution is outside
this skill and needs an explicit device-specific request and operational review.

See the official [visualizer example](https://docs.pylabrobot.org/dev/user_guide/machine-agnostic-features/using-the-visualizer.html)
for chatterbox setup; use documentation matching the selected release when
running it. No robot connection is needed to consult or apply this skill.

Migration source: Biomni commit `400c1f366b96a35ca253e13c9b06c5076af41d65`,
`biomni/tool/lab_automation.py`. Its source-rewriting, import-execution and
thread-timeout testing wrapper is intentionally not retained.

Run the reviewed simulation script with `biov exec python -- simulation.py`.
Keep the explicit simulation backend and inspect its actual output.

Read and write files in the working directory or configured data directories.
