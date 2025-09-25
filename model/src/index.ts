import type { GraphMakerState } from '@milaboratories/graph-maker';
import type {
  InferOutputsType,
  PFrameHandle,
  PlDataTableStateV2,
  PlRef,
} from '@platforma-sdk/model';
import {
  BlockModel,
  createPFrameForGraphs,
  createPlDataTableSheet,
  createPlDataTableStateV2,
  createPlDataTableV2,
  getUniquePartitionKeys,
  isPColumnSpec,
} from '@platforma-sdk/model';

export type UiState = {
  graphStateBubble: GraphMakerState;
  graphStateUMAP: GraphMakerState;
  graphStateTSNE: GraphMakerState;
  tableState: PlDataTableStateV2;
};

export type BlockArgs = {
  clusterAnnotationRef?: PlRef;
  title?: string;
  topN: number;
  logfcCutoff: number;
  pvalCutoff: number;
  strictOverlap: boolean;
};

export const model = BlockModel.create()

  .withArgs<BlockArgs>({
    topN: 3,
    logfcCutoff: 1.0,
    pvalCutoff: 0.01,
    strictOverlap: false,
  })

  .argsValid((ctx) => {
    // Check if cluster annotation is selected
    if (!ctx.args.clusterAnnotationRef) {
      return false;
    }

    // Check if topN is a valid positive number
    if (!ctx.args.topN || ctx.args.topN < 1) {
      return false;
    }

    return true;
  })

  .withUiState<UiState>({
    graphStateBubble: {
      title: 'Dotplot',
      template: 'bubble',
      layersSettings: {
        bubble: {
          normalizationDirection: null,
        },
      },
    },
    graphStateUMAP: {
      title: 'UMAP',
      template: 'dots',
    },
    graphStateTSNE: {
      title: 'tSNE',
      template: 'dots',
    },
    tableState: createPlDataTableStateV2(),
  })

  // Allow inputs from any single-cell grouping block
  .output('clusterAnnotationOptions', (ctx) =>
    ctx.resultPool.getOptions((spec) => isPColumnSpec(spec)
      && (spec.name === 'pl7.app/rna-seq/leidencluster'
        || spec.name === 'pl7.app/rna-seq/cellType')
    , { includeNativeLabel: true, addLabelAsSuffix: true }),
  )

  .output('clusterMarkersPt', (ctx) => {
    const pCols = ctx.outputs?.resolve('clusterMarkersPf')?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }

    return createPlDataTableV2(ctx, pCols, ctx.uiState.tableState);
  })

  .output('clusterMarkersSheets', (ctx) => {
    const pCols = ctx.outputs?.resolve('clusterMarkersPf')?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }

    const anchor = pCols[0];
    if (!anchor) return undefined;

    const r = getUniquePartitionKeys(anchor.data);
    if (!r) return undefined;

    return r.map((values, i) => createPlDataTableSheet(ctx, anchor.spec.axesSpec[i], values));
  })

  .output('clusterMarkersTopPf', (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve('clusterMarkersTopPf')?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }
    return createPFrameForGraphs(ctx, pCols);
  })

  .output('umapPf', (ctx): PFrameHandle | undefined => {
    return createPFrameForGraphs(ctx);
  })

  .output('isRunning', (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .sections((_ctx) => ([
    { type: 'link', href: '/', label: 'Main' },
    // { type: 'link', href: '/umap', label: 'UMAP' },
    { type: 'link', href: '/dotplot', label: 'Dotplot' },
  ]))

  .title((ctx) =>
    ctx.args.title
      ? `Cluster Markers - ${ctx.args.title}`
      : 'Cluster Markers',
  )

  .done();

export type BlockOutputs = InferOutputsType<typeof model>;
