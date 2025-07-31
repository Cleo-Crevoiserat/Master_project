*Forall-Exist Statements & Fair Allocation problems*

This project serves **two purposes**:

1. **Logical Verification**: Given a matrix \( W \in \mathbb{Z}^{m 	imes n} \) and a polyhedron \( Q \subseteq \mathbb{R}^m \), verify whether:
   \[
   orall b \in Q \cap \mathbb{Z}^m, \ \exists x \in \mathbb{Z}^n 	ext{ such that } Wx \leq b
   \]

2. **Fair Allocation**: For \( n \) agents and \( m \) object types, determine whether a **fair allocation** exists, and if so, compute it.

---

## ⚙️ Requirements

- **Gurobi** optimizer (make sure it's installed and accessible with a valid license)
- C++ compiler (e.g., `g++`, `clang++`)

---

## 🧪 Usage – Case 1: Forall–Exist Verification

In your `main()`:
1. Set `prob = 0`
2. Define matrix dimensions `m`, `n`
3. Initialize matrix \( W \) using Eigen:

```cpp
W << 1, 2,
     3, 6,
     5, 7;
```

4. Define polyhedron \( Q \) as a list of inequalities in the form:

```
Q = {
  {c, a₁, a₂, ..., aₙ},
  {d, b₁, b₂, ..., bₙ},
  ...
};
```

Each row corresponds to the inequality:

\[
a_1 x_1 + a_2 x_2 + \dots + a_n x_n + c \geq 0
\]

---

## 🎯 Usage – Case 2: Fair Allocation

Set up the problem with:
- `m`: number of object types
- `n`: number of agents
- `N`: vector indicating how many objects of each type (e.g., `N << 1, 2, 3, 4, 5;`)
- Utility function for each agent: a `vector<int>` of length `m` (e.g., `{3, 2, 0, 1, 4}`)

Each agent's utility is defined per object type.

---

## 📎 Notes

- The code leverages Gurobi for solving ILPs; make sure your environment is correctly configured.
- The logic is inspired by recent research (e.g., on pseudopolynomial forall-exist algorithms).
