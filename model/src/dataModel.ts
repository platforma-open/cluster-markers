import { kind } from "@platforma-open/milaboratories.cluster-markers.kind";
import { createPlDataTableStateV2, DataModelBuilder } from "@platforma-sdk/model";
import type { BlockArgs, BlockData, BlockUiState } from "./types";

export const blockDataModel = new DataModelBuilder({ kind })
  .from<BlockData>("v1")
  // V1 kept these in two buckets and V3 keeps one, but the fields are the same
  // fields — the merge is the whole upgrade.
  .upgradeLegacy<BlockArgs, BlockUiState>(({ args, uiState }) => ({
    ...args,
    ...uiState,
    // Not on disk under V1 — the sheet was the only way to read the table, so
    // that is what a migrated project keeps seeing.
    tableScope: "cluster" as const,
  }))
  .init(({ params }) => ({
    clusterAnnotationRef: params?.clusterAnnotationRef,
    topN: params?.topN ?? 3,
    logfcCutoff: params?.logfcCutoff ?? 1.0,
    pvalCutoff: params?.pvalCutoff ?? 0.01,
    strictOverlap: params?.strictOverlap ?? false,

    tableScope: "cluster" as const,

    graphStateBubble: {
      title: "Dotplot",
      template: "bubble",
      layersSettings: {
        bubble: {
          normalizationDirection: null,
        },
      },
    },
    graphStateUMAP: {
      title: "UMAP",
      template: "dots",
    },
    graphStateTSNE: {
      title: "tSNE",
      template: "dots",
    },
    tableState: createPlDataTableStateV2(),
  }));
