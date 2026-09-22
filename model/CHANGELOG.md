# @platforma-open/milaboratories.cluster-markers.model

## 1.6.0

### Minor Changes

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

- d6bce9a: Add a per-cluster / all-clusters scope switch to the marker table.

  The table's Export button writes the table you are looking at, and the cluster
  axis was a sheet — so exporting gave one cluster's markers, without a cluster
  column. The switch beside the cluster picker drops the sheet, putting every
  cluster in the table at once with Cluster as an ordinary column, so Export
  writes the whole result set.

  Each scope keeps its own sorting, filters and column layout. The default is
  per-cluster, which is what the block did before.

### Patch Changes

- Updated dependencies [3319e63]
  - @platforma-open/milaboratories.cluster-markers.kind@1.0.1

## 1.5.0

### Minor Changes

- d2368a9: Enable block deduplication and update metadata.

## 1.4.4

### Patch Changes

- 8b7a188: technical release
- 99978bc: technical release
- 84bb66a: technical release
- 89db503: technical release

## 1.4.3

### Patch Changes

- bae60fb: Support any cell group, not only leiden clusters

## 1.4.2

### Patch Changes

- 4c6b675: Update SDK

## 1.4.1

### Patch Changes

- 75d2a9c: Fixed github build

## 1.4.0

### Minor Changes

- 42f571e: Update SDK packages, minor plot fixes, migrate to PlAgDataTableV2 and expose option for cluster marker overlap filtering

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
