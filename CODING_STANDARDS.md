# MRCPP Coding Standards

Conventions that Claude (and contributors) must follow when writing **new** code in this repository.

## Scope

- Applies to new files and to new or modified lines. Do **not** reformat or rename legacy code as part of unrelated changes; legacy code will be migrated gradually in dedicated PRs.
- When a rule cannot be applied locally without making a legacy file inconsistent (for example adding a method to a class whose other members use a different scheme), match the surrounding code and leave the migration for later.
- Renaming anything in the public API (`api/`, public headers) is a breaking change for downstream users such as MRChem. Do not do it without an explicit decision.
- Background: see `docs/CODING_STANDARDS_AUDIT.md` for the survey of the existing conventions.

## Naming

| Entity | Convention | Example |
|---|---|---|
| Class names | `PascalCase` | `FunctionTree` |
| Class file names | `PascalCase` | `FunctionTree.h`, `FunctionTree.cpp` |
| Member functions | `camelCase` | `getScale()`, `isLeafNode()` |
| Non-member functions | `snake_case` | `apply_far_field()`, `compute_norm()` |
| Member variables | `camelCase` | `squareNorm`, `maxDepth` |
| Global variables | `SNAKE_CASE` (upper case) | `MACHINE_ZERO`, `MAX_DEPTH` |
| Local variables | `snake_case` | `node_index`, `max_scale` |

Notes:

- "Non-member functions" covers free functions and functions in utility namespaces (`mrcpp::math_utils`, `mrcpp::tree_utils`, ...).
- `SNAKE_CASE` for globals means upper-case words separated by underscores.
- Files that hold only non-member functions (not a class) are named in `snake_case`, as the existing `project.h`, `apply.h`, `math_utils.h` already are.

## Control-flow style (braces and indentation)

Applies to `if` / `else if` / `else`, `for`, `while`, `do ... while` and similar blocks.

Three forms are allowed:

```cpp
if (a) do_b();      // allowed (single short statement on the same line)

if (a)              // allowed, preferred
    do_b();

if (a) {            // allowed, preferred
    do_b();
}
```

- The two **preferred** forms are the indented single statement without braces and the braced block with the opening brace attached to the statement (`) {`).
- Use the braced form when the body has more than one statement, or is long or complex.
- The one-line form is for trivial statements only.
- Indent with 4 spaces, never tabs.
- The same three forms apply to the other blocks:

```cpp
for (auto i = 0; i < n; i++) sum += x[i];   // allowed

for (auto i = 0; i < n; i++)                // allowed, preferred
    sum += x[i];

for (auto i = 0; i < n; i++) {              // allowed, preferred
    sum += x[i];
}

if (a) {
    do_b();
} else {
    do_c();
}
```
