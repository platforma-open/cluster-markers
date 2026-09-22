import type { GraphMakerState } from "@milaboratories/graph-maker";
import type { PlDataTableStateV2, PlRef } from "@platforma-sdk/model";

/**
 * Workflow-facing args — every user-committed analysis decision. Unchanged
 * from the V1 shape on purpose: the workflow reads these keys verbatim, so
 * keeping them byte-identical is what stops every existing project going
 * stale the moment it is opened on this version.
 *
 * `title` is never written by the UI; it survives here so `.title()` keeps
 * reading the same field it always did.
 */
export type BlockArgs = {
  clusterAnnotationRef?: PlRef;
  title?: string;
  topN: number;
  logfcCutoff: number;
  pvalCutoff: number;
  strictOverlap: boolean;
};

/**
 * Which slice of the marker table the Main page shows.
 */
export type TableScope = "cluster" | "all";

/** The view state the workflow never sees. Read only by `.upgradeLegacy` under this name. */
export type BlockUiState = {
  graphStateBubble: GraphMakerState;
  graphStateUMAP: GraphMakerState;
  graphStateTSNE: GraphMakerState;
  tableState: PlDataTableStateV2;
};

/**
 * Unified V3 data — the args above plus the view state: the four plot/table
 * states and the scope the marker table is read at.
 *
 * The intersection is the honest shape here, not V1 in a V3 jacket. Every
 * settings-panel field is an analysis decision the user commits and then
 * presses Run for, so none of them belongs in `prerunArgs`; this block
 * discovers nothing and stages nothing. And no V1 field was bent to dodge the
 * staleness gate, so there is no distortion to undo — the split that already
 * existed is the split V3 wants.
 */
export type BlockData = BlockArgs &
  BlockUiState & {
    tableScope: TableScope;
  };
