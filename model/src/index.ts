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
    // I've added these "||" for backward compatibility (As I see, the shape of PColum was changed)
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
    // let pCols = ctx.outputs?.resolve('clusterMarkersPf')?.getPColumns();
    const pCols = ctx.outputs?.resolve('clusterMarkersPf')?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }

    // Filter by selected comparison
    // pCols = pCols.filter(
    //   (col) => col.spec.axesSpec[0]?.domain?.['pl7.app/comparison'] === ctx.uiState.comparison,
    // );

    // return createPlDataTable(ctx, pCols, ctx.uiState?.tableState);

    const anchor = pCols[0];
    if (!anchor) return undefined;

    const r = getUniquePartitionKeys(anchor.data);
    if (!r) return undefined;

    return {
      table: createPlDataTable(ctx, pCols, ctx.uiState?.tableState),
      sheets: r.map((values, i) => createPlDataTableSheet(ctx, anchor.spec.axesSpec[i], values)),
    };
  })

// .output('clusterMarkersPf', (ctx): PFrameHandle | undefined => {
//   const pCols = ctx.outputs?.resolve('clusterMarkersPf')?.getPColumns();
//   if (pCols === undefined) {
//     return undefined;
//   }
// })

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

// .output('UMAPPf', (ctx): PFrameHandle | undefined => {
//   const pCols = ctx.outputs?.resolve('UMAPPf')?.getPColumns();
//   if (pCols === undefined) {
//     return undefined;
//   }

//   // enriching with upstream data
//   const upstream = ctx.resultPool
//     .getData()
//     .entries.map((v) => v.obj)
//     .filter(isPColumn)
//     .filter((column) => column.id.includes('metadata'));

//   return ctx.createPFrame([...pCols, ...upstream]);
// })

// .output('tSNEPf', (ctx): PFrameHandle | undefined => {
//   const pCols = ctx.outputs?.resolve('tSNEPf')?.getPColumns();
//   if (pCols === undefined) {
//     return undefined;
//   }

//   // enriching with upstream data
//   const upstream = ctx.resultPool
//     .getData()
//     .entries.map((v) => v.obj)
//     .filter(isPColumn)
//     .filter((column) => column.id.includes('metadata'));

//   return ctx.createPFrame([...pCols, ...upstream]);
// })

  .sections((_ctx) => ([
    { type: 'link', href: '/', label: 'Main' },
    { type: 'link', href: '/umap', label: 'UMAP' },
    { type: 'link', href: '/dotplot', label: 'Dotplot' },
  ]))

  .done();

export type BlockOutputs = InferOutputsType<typeof model>;
