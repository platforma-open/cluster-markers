import type { GraphMakerState } from '@milaboratories/graph-maker';
import type {
  InferOutputsType,
  PFrameHandle,
  PlDataTableState,
  PlRef } from '@platforma-sdk/model';
import {
  BlockModel,
  createPFrameForGraphs,
  createPlDataTable,
  createPlDataTableSheet,
  getUniquePartitionKeys,
  isPColumn,
  // isPColumn,
  isPColumnSpec,
} from '@platforma-sdk/model';

export type UiState = {
  graphStateBubble: GraphMakerState;
  graphStateUMAP: GraphMakerState;
  graphStateTSNE: GraphMakerState;
  tableState: PlDataTableState;
};

export type BlockArgs = {
  countsRef?: PlRef;
  clusterAnnotationRef?: PlRef;
  title?: string;
};

export const model = BlockModel.create()

  .withArgs<BlockArgs>({
  })

  .withUiState<UiState>({
    graphStateBubble: {
      title: 'Dotplot',
      template: 'bubble',
    },
    graphStateUMAP: {
      title: 'UMAP',
      template: 'dots',
    },
    graphStateTSNE: {
      title: 'tSNE',
      template: 'dots',
    },
    tableState: {
      gridState: {},
      pTableParams: {
        sorting: [],
        filters: [],
      },
    },
  })

  .output('countsOptions', (ctx) =>
    ctx.resultPool.getOptions((spec) => isPColumnSpec(spec)
      && spec.name === 'pl7.app/rna-seq/countMatrix' && spec.domain?.['pl7.app/rna-seq/normalized'] === 'false'
    , { includeNativeLabel: true, addLabelAsSuffix: true }),
  )

  .output('clusterAnnotationOptions', (ctx) =>
    ctx.resultPool.getOptions((spec) => isPColumnSpec(spec)
      && spec.name === 'pl7.app/rna-seq/leidencluster'
    , { includeNativeLabel: true, addLabelAsSuffix: true }),
  )

  .output('clusterMarkersPt', (ctx) => {
    const pCols = ctx.outputs?.resolve('clusterMarkersPf')?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }

    const anchor = pCols[0];
    if (!anchor) return undefined;

    const r = getUniquePartitionKeys(anchor.data);
    if (!r) return undefined;

    return {
      table: createPlDataTable(ctx, pCols, ctx.uiState?.tableState),
      sheets: r.map((values, i) => createPlDataTableSheet(ctx, anchor.spec.axesSpec[i], values)),
    };
  })

  .output('clusterMarkersTop3Pf', (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve('clusterMarkersTop3Pf')?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }
    return createPFrameForGraphs(ctx, pCols);
  })

  .output('UMAPPf', (ctx): PFrameHandle | undefined => {
    return createPFrameForGraphs(ctx,
      ctx.resultPool
        .getData()
        .entries.map((c) => c.obj)
        .filter(isPColumn)
        .filter((column) => column.spec.name.includes('umap')),
    );
  })

  .output('isRunning', (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .sections((_ctx) => ([
    { type: 'link', href: '/', label: 'Main' },
    { type: 'link', href: '/umap', label: 'UMAP' },
    { type: 'link', href: '/dotplot', label: 'Dotplot' },
  ]))

  .title((ctx) =>
    ctx.args.title
      ? `Cluster Markers - ${ctx.args.title}`
      : 'Cluster Markers',
  )

  .done();

export type BlockOutputs = InferOutputsType<typeof model>;
