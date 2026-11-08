title: Related Pages
summary: Supplementary documentation, theoretical guides, and integration tutorials for MUSCLE.

# Supplementary Documentation

Welcome to the extended documentation for **MUSCLE**. This section provides in-depth information that complements the technical API reference generated from the source code.

Select a topic below to learn more:

### [Theoretical Background](theory.html)
Mathematical foundation of the library, including plastic return mapping, yield criteria (Von Mises, CPB06), and hardening laws.

### [Tensor Algebra Guide](tensors.html)
A guide on using the optimized tensor engine, custom operators (`.ddot.`, `.tdotsym.`), and Voigt notation storage.

### [FEA Integration](integration.html)
How to call a declared material from the user subroutines of Abaqus (UMAT), ANSYS (USERMAT) and LS-DYNA (umat41).

### [Developer Guide](contributing.html)
How to add a model, a material or a solver, with checklists.

### [Conventions](conventions.html)
The rules every contribution follows: code style, `sp`/`dp`, 2D/3D, minimal FYPP, GPU-safe code and tests.

### [Architecture](architecture.html)
How a declared material is generated from the registry, the solver bodies and the model families.

### [Glossary](glossary.html)
Abbreviations of the names of materials, models and solvers (`cp`, `pz`, `vm`, `hsv`, ...).

### [Migration from the OO API](migration.html)
Equivalence table from the former object-oriented API (`Elasticity_linear`, `Closest_point`, `Plastic_material_history`, ...) to the declared materials, and what is no longer available.