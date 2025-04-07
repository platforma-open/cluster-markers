<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlBlockPage } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';

import { GraphMaker } from '@milaboratories/graph-maker';
import type { PredefinedGraphOption } from '@milaboratories/graph-maker';
import { watch } from 'vue';

const app = useApp();

watch(() => app.model.ui, async (value) => {
  console.log(value, 'app.model.ui.graphStateBubble');
}, { immediate: true });

const defaultOptions: PredefinedGraphOption<'bubble'>[] = [
  {
    inputName: 'valueSize',
    selectedSource: {
      kind: 'PColumn',
      name: 'pl7.app/rna-seq/percentcells',
      valueType: 'Double',
      axesSpec: [
        {
          name: 'pl7.app/rna-seq/cluster-num',
          type: 'String',
        },
        {
          name: 'pl7.app/rna-seq/geneId',
          type: 'String',
        },
      ],
    },
  },
  {
    inputName: 'valueColor',
    selectedSource: {
      kind: 'PColumn',
      name: 'pl7.app/rna-seq/meanexpression',
      valueType: 'Double',
      axesSpec: [
        {
          name: 'pl7.app/rna-seq/cluster-num',
          type: 'String',
        },
        {
          name: 'pl7.app/rna-seq/geneId',
          type: 'String',
        },
      ],
    },
  },
  {
    inputName: 'x',
    selectedSource: {
      name: 'pl7.app/rna-seq/cluster-num',
      type: 'String',
    },
  },
  {
    inputName: 'y',
    selectedSource: {
      name: 'pl7.app/rna-seq/geneId',
      type: 'String',
    },
  },
];

</script>

<template>
  <PlBlockPage>
    <GraphMaker
      v-model="app.model.ui.graphStateBubble" chartType="bubble"
      :p-frame="app.model.outputs.clusterMarkersTop3Pf" :default-options="defaultOptions"
    />
  </PlBlockPage>
</template>
