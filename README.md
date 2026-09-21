# This is my quantum mechanics course

Quantum Mechanics Problems: A Computational Approach
https://img.shields.io/badge/Made%20with-Python-1f425f.svg
https://img.shields.io/badge/University-[Your%20University%20Name]-blue
https://img.shields.io/badge/Course-Quantum%20Mechanics-orange

📖 Overview
This repository contains a collection of Python solutions and simulations for various problems encountered in a university-level Quantum Mechanics course.

The goal of this project is to bridge the gap between abstract theoretical physics and computational implementation. By using Python, complex quantum phenomena—which are often difficult to visualize analytically—are simulated, visualized, and analyzed numerically.

🎓 Course Information
Course Name: [Exercises of Quantum Mechanics]

Institution: [Aristotle University of Thessaloniki]

Semester: [Fall 2026]

Instructor: [Prof. Diakonidis Theodoros]

🛠️ Technologies Used
The simulations and solutions in this repository rely on the following Python libraries:

NumPy: For efficient matrix operations and numerical calculations.

SciPy: For solving differential equations (Schrödinger equation) and linear algebra routines.

Matplotlib: For visualizing wavefunctions, probability densities, and energy levels.

SymPy (Optional): Used for symbolic derivation of quantum operators.

📂 Repository Structure
The problems are organized by topic or chapter. Each directory contains the Python script/Jupyter Notebook, the resulting plots, and a brief explanation of the physics involved.

├── 01_Wave_Mechanics/
│   ├── particle_in_a_box.py       # 1D and 2D Infinite Square Well
│   ├── finite_square_well.py      # Bound states in a finite potential
│   └── quantum_harmonic.py        # Harmonic Oscillator wavefunctions
├── 02_Formalism/
│   ├── dirac_notation.py          # Bra-ket manipulation
│   └── commutators.py             # Verifying uncertainty principles
├── 03_Angular_Momentum/
│   ├── spherical_harmonics.py     # Visualizing Y_lm
│   └── spin_1_2.py               # Pauli matrices and spin states
├── 04_Approximation_Methods/
│   ├── variational_principle.py   # Estimating ground state energy
│   └── wkb_approximation.py       # Tunneling probabilities
└── README.md
