# @platforma-open/milaboratories.cluster-markers.kind

## 1.0.1

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
