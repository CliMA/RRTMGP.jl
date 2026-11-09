# Solver

```@meta
CurrentModule = RRTMGP.RTE
```

```@docs
NoScatLWRTE
TwoStreamLWRTE
NoScatSWRTE
TwoStreamSWRTE
```

## G-point groups

A two-stream solve can spread its g-point loop across the threads of one wider
kernel instead of running it on one thread per column. These hold the per-group
workspaces that make that possible.

```@docs
LWGptGroups
SWGptGroups
```
