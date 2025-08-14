# Forall-Exist Statements & Fair Allocation problems

This project has **two different purposes**:

1. **Solve a forall exist statement**: Given a matrix $( W \in \mathbb{Z}^{m 	\times n} )$ and a polyhedron $( Q \subseteq \mathbb{R}^m )$, decide the validity of the statement:
   $$\forall b \in Q \cap \mathbb{Z}^m,  \exists x \in \mathbb{Z}^n 	\text{ such that } Wx \leq b$$

2. **Solve a Fair Allocation problem**: For $n$ agents and $m$ object categories, determine whether a **fair allocation** exists, and if so, compute it.

---

## Requirements

- **Gurobi** optimizer (make sure it's installed and accessible with a valid license)
- C++ compiler (e.g., `g++`, `clang++`)

---
# How to use it
 Case 1: Forall–Exist Verification

In the `main()`:
1. Set `prob = 0`
2. Define matrix dimensions `m`, `n`
3. Initialize matrix $W$ using Eigen.  
 Example:  
  ![Matrice exemple](docs/images/matrice_exemple.png)
```cpp
W << 1, 2,
     3, 6,
     5, 7;
```

4. Define polyhedron  $Q$ as a list of inequalities in the form:

```
Q = {
  {c, a₁, a₂, ..., aₙ},
  {d, b₁, b₂, ..., bₙ},
  ...
};
```

Each row corresponds to the inequality:

$a_1 x_1 + a_2 x_2 + \dots + a_n x_n + c \geq 0$

---

## Case 2: Fair Allocation

Set up the problem with:
- `m`: number of object types
- `n`: number of agents
- `N`: vector indicating how many objects of each category (e.g., `N << 1, 2, 3, 4, 5;` if you have one item for the first category, two for the second one and so on.)
- Utility function for each agent: a `vector<int>` of length `m` (e.g., `{3, 2, 0, 1, 4}` means that $U_{a}(1)=3$, $U_{a}(2)=2$...)

Each agent's utility is defined per object type.

---

## 📎 Notes

- The code leverages Gurobi for solving ILPs; make sure your environment is correctly configured.

## Structure
### Code
Everything related to code is in the `src` folder. This includes `.cpp` files and their headers as well as the compiled exectutable scripts.

The script `Forall-exist.exe` enables to deploy the tools we developed in the project to solve both a fair allocation and forall-exists statement.

### Docs
The images used in the readme are in the `images` folder.


## Running code
If you want to directy run the program, you can simply use the compiled `.exe` file. Navigate to the `src` folder and run the following command:

```sh
Forall-exist.exe <problem_type> <size> 2  1 1  0 0 1 1
```
## Compiling C++ code
If you need to modify the code or re-compile it for some reason, you can navigate to the `src` folder and run something along the line.

```sh
cl /EHsc /MD /std:c++17 /Iheaders /I"C:\gurobi1103\win64\include" /I"C:\Libraries\eigen-3.4.0" cpp\Forall_Exist.cpp cpp\Construction_Allocation.cpp cpp\Construction_Part_1.cpp cpp\Solve_6.cpp cpp\Useful_fct.cpp C:\gurobi1103\win64\lib\gurobi110.lib C:\gurobi1103\win64\lib\gurobi_c++md2017.lib /FeForall-exist.exe && Forall-exist.exe
```

*Note: On windows, you might want to run the command on a `x64 Native Tools Command Prompt` to be able to run `cl`.*