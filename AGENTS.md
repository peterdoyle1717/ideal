# Repo rules for Codex

## Speed discipline (HARD)

Production solver sweeps run in C, link sparse SuperLU, and process a
chunk of CLERS per worker process. Do not spawn one process per CLERS
or one Python interpreter per CLERS/net/case.

Python may orchestrate coarse chunks, prepare inputs, or run one-off
experiments, but it must not appear in a sweep's hot loop. If you are
writing `subprocess.run([python, ...])`, `python script.py`, or
equivalent inside a per-CLERS/per-net loop, stop and move that loop
into a C binary.

Before launching any sweep, state:

- the binary that will run the hot loop
- where the sparse linear solve happens
- the chunk size
- the total number of worker processes spawned for N CLERS

If spawned processes are approximately N, or the hot-loop binary is a
Python interpreter, that is the known failure mode and the sweep must
not be launched.

Shell wrappers count as the binary they ultimately invoke; inspect
them before launch.

Exception: calling a fast C/C++ helper once per item is allowed only
when that helper is not the solver bottleneck and does not itself
launch Python. Name the helper and justify why it is outside the hot
solver path.

## Reviewer behaviour

When reviewing a sweep harness, check the speed-discipline rule
mechanically:

- Does the per-CLERS code path call a Python interpreter? If yes,
  REQUEST CHANGES.
- How many processes are spawned for N CLERS? If ≈ N, REQUEST CHANGES.
- Where does the sparse linear solve happen? It must be in a C binary
  linked against SuperLU.
