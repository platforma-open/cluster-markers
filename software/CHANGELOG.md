# @platforma-open/milaboratories.cluster-markers.software

## 1.6.7

### Patch Changes

- 3319e63: Adopt the block-tools structure, add the block kind, and move the model to
  BlockModelV3.

  The block gains an init-params contract — the cluster annotation, the
  top-markers count, the two significance thresholds and the specificity mode —
  so a project template can create it pre-configured, and `platforma` replaces
  `model` as the model's export. Persisted settings carry over unchanged through
  the legacy upgrader, and the args the workflow reads are the same keys as
  before, so existing projects do not go stale.

  Also a full SDK upgrade: model 1.49 to 1.83, workflow-tengo 5.7 to 6.11 and
  tengo-builder 2.4 to 4.1.

## 1.6.6

### Patch Changes

- 8ee5a70: Improve performance

## 1.6.5

### Patch Changes

- 8b7a188: technical release
- 99978bc: technical release
- 84bb66a: technical release
- 89db503: technical release

## 1.6.4

### Patch Changes

- bae60fb: Support any cell group, not only leiden clusters

## 1.6.3

### Patch Changes

- fae48db: Full SDK update

## 1.6.2

### Patch Changes

- 75d2a9c: Fixed github build

## 1.6.1

### Patch Changes

- b0213cd: Updated SDK and added running bar

## 1.6.0

### Minor Changes

- 5ca4a45: Update script to work with updated Cell ID axis specifications

## 1.5.0

### Minor Changes

- 72dc2fe: allow create venv on Windows

## 1.4.0

### Minor Changes

- 4c10779: Add DEG list exports for functional analysis and filtering

## 1.3.0

### Minor Changes

- 5d13c7b: Included arguments to select Fold Change and P-Value cutoffs

## 1.2.0

### Minor Changes

- fc0481c: Added option to select number of marker genes to be used in visualizations

## 1.1.1

### Patch Changes

- 09cb764: Improved visualization

## 1.1.0

### Minor Changes

- c25661f: First version
