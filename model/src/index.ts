import { kind } from "@platforma-open/milaboratories.cluster-markers.kind";
import type { InferOutputsType, PFrameHandle } from "@platforma-sdk/model";
import {
  BlockModelV3,
  createPFrameForGraphs,
  createPlDataTableSheet,
  createPlDataTableV2,
  getUniquePartitionKeys,
  isPColumnSpec,
} from "@platforma-sdk/model";
import { blockDataModel } from "./dataModel";
import type { BlockArgs } from "./types";

export { blockDataModel } from "./dataModel";
export type { BlockArgs, BlockData, BlockUiState } from "./types";

export const platforma = BlockModelV3.create({ dataModel: blockDataModel, kind })

  // The run gate. Throwing surfaces the reason in the UI, which a disabled Run
  // button would not. Replaces V1's `.argsValid`, so validation cannot drift
  // from the projection any more.
  .args<BlockArgs>((data): BlockArgs => {
    const {
      graphStateBubble: _graphStateBubble,
      graphStateUMAP: _graphStateUMAP,
      graphStateTSNE: _graphStateTSNE,
      tableState: _tableState,
      tableScope: _tableScope,
      ...args
    } = data;

    if (args.clusterAnnotationRef === undefined) {
      throw new Error("Cluster annotation is required");
    }
    if (!args.topN || args.topN < 1) {
      throw new Error("Top markers per cluster must be at least 1");
    }

    // Nothing is canonicalized here on purpose: these keys are what the V1
    // workflow already read, and the workflow still supplies its own defaults
    // for the two cutoffs a user can clear. Substituting them here would
    // change the args bytes of every project that has one cleared, and stale
    // it for no gain.
    return args;
  })

  // Inverse of the kind's init-params contract: the five settings a user picks
  // by hand. `title` is left out — nothing writes it, so there is nothing to
  // carry into a template.
  .templateParams((data) => ({
    clusterAnnotationRef: data.clusterAnnotationRef,
    topN: data.topN,
    logfcCutoff: data.logfcCutoff,
    pvalCutoff: data.pvalCutoff,
    strictOverlap: data.strictOverlap,
  }))

  // Allow inputs from any single-cell grouping block
  .output("clusterAnnotationOptions", (ctx) =>
    ctx.resultPool.getOptions(
      (spec) =>
        isPColumnSpec(spec) &&
        (spec.name === "pl7.app/rna-seq/leidencluster" || spec.name === "pl7.app/rna-seq/cellType"),
      { includeNativeLabel: true, addLabelAsSuffix: true },
    ),
  )

  // `withStatus` on the next three: PlAgDataTableV2 (via usePlDataTableSettingsV2)
  // and GraphMaker both take the status-wrapped output, and render the pending
  // and error states from it themselves.
  .outputWithStatus("clusterMarkersPt", (ctx) => {
    const pCols = ctx.outputs?.resolve("clusterMarkersPf")?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }

    return createPlDataTableV2(ctx, pCols, ctx.data.tableState);
  })

  .output("clusterMarkersSheets", (ctx) => {
    const pCols = ctx.outputs?.resolve("clusterMarkersPf")?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }

    const anchor = pCols[0];
    if (!anchor) return undefined;

    const r = getUniquePartitionKeys(anchor.data);
    if (!r) return undefined;

    return r.map((values, i) => createPlDataTableSheet(ctx, anchor.spec.axesSpec[i], values));
  })

  .outputWithStatus("clusterMarkersTopPf", (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve("clusterMarkersTopPf")?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }
    return createPFrameForGraphs(ctx, pCols);
  })

  .outputWithStatus("umapPf", (ctx): PFrameHandle | undefined => {
    return createPFrameForGraphs(ctx);
  })

  .output("isRunning", (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .sections((_ctx) => [
    { type: "link", href: "/", label: "Main" },
    // { type: 'link', href: '/umap', label: 'UMAP' },
    { type: "link", href: "/dotplot", label: "Dotplot" },
  ])

  .title((ctx) => (ctx.data.title ? `Cluster Markers - ${ctx.data.title}` : "Cluster Markers"))

  .done();

export type BlockOutputs = InferOutputsType<typeof platforma>;
