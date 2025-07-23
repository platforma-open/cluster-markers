<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAgDataTableV2, PlAlert, PlBlockPage, PlBtnGhost, PlBtnGroup, PlDropdownRef, PlMaskIcon24, PlNumberField, PlRow, PlSlideModal, usePlDataTableSettingsV2 } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';
import { reactive } from 'vue';

const app = useApp();

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.clusterMarkersPt,
  sheets: () => app.model.outputs.clusterMarkersSheets,
});

const overlapOptions = [
  { text: 'Non-exclusive', value: false },
  { text: 'Strict overlap', value: true },
];

const data = reactive<{
  settingsOpen: boolean;
}>({
  settingsOpen: app.model.args.clusterAnnotationRef === undefined,
});

</script>

<template>
  <PlBlockPage>
    <template #title>Cluster Marker Discovery</template>
    <template #append>
      <PlBtnGhost @click.stop="() => data.settingsOpen = true">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>
    <PlAgDataTableV2
      v-model="app.model.ui.tableState"
      :settings="tableSettings"
      show-export-button
    />
    <PlSlideModal v-model="data.settingsOpen">
      <template #title>Settings</template>
      <PlDropdownRef
        v-model="app.model.args.clusterAnnotationRef" :options="app.model.outputs.clusterAnnotationOptions"
        label="Cluster annotation"
        clearable
      />
      <PlNumberField
        v-model="app.model.args.topN"
        label="Top markers per cluster" :minValue="1" :step="1"
      >
        <template #tooltip>
          <div>
            <strong>Number of top markers to display</strong><br/>
            Determines how many of the most significant marker genes will be shown for each cluster in the results table. Higher values provide more comprehensive marker profiles but may include less specific markers.
          </div>
        </template>
      </PlNumberField>
      <PlBtnGroup
        v-model="app.model.args.strictOverlap"
        label="Marker specificity mode"
        :options="overlapOptions"
      >
        <template #tooltip>
          <div>
            <strong>Marker specificity filtering</strong><br/>
            Controls how strictly markers are filtered based on their expression specificity to clusters.<br/><br/>
            <strong>Non-exclusive:</strong> Includes genes expressed in ≥20% of cells within the target cluster, regardless of expression in other clusters. More permissive, captures broader marker profiles.<br/><br/>
            <strong>Strict overlap:</strong> Requires genes to be expressed in ≥20% of cells in the target cluster AND &lt;20% of cells in other clusters. More stringent, identifies highly cluster-specific markers.<br/><br/>
          </div>
        </template>
      </PlBtnGroup>
      <PlRow>
        <PlNumberField
          v-model="app.model.args.logfcCutoff"
          label="Log2(FC)" :minValue="0" :step="0.1"
        >
          <template #tooltip>
            <div>
              <strong>Log2 fold-change threshold</strong><br/>
              Minimum absolute log2 fold-change required for a gene to be considered a cluster marker. Higher values (e.g., 1.0-2.0) identify more dramatically upregulated genes, while lower values (e.g., 0.5-0.8) capture subtler expression differences.
            </div>
          </template>
        </PlNumberField>
        <PlNumberField
          v-model="app.model.args.pvalCutoff"
          label="Adjusted p-value" :minValue="0" :maxValue="1" :step="0.01"
        >
          <template #tooltip>
            <div>
              <strong>Statistical significance threshold</strong><br/>
              Maximum adjusted p-value for marker gene significance. Standard values are 0.05 (5% false discovery rate) or 0.01 (1% false discovery rate) for more stringent filtering.
            </div>
          </template>
        </PlNumberField>
      </PlRow>
      <!-- Add warnings if selected threshold are out of most commonly used bounds -->
      <PlAlert v-if="app.model.args.pvalCutoff > 0.05" type="warn">
        {{ "Warning: The selected adjusted p-value threshold is higher than the commonly recommended 0.05" }}
      </PlAlert>
      <PlAlert v-if="app.model.args.logfcCutoff < 0.6" type="warn">
        {{ "Warning: The selected Log2(FC) threshold may be too low for identifying robust cluster markers" }}
      </PlAlert>
    </PlSlideModal>
  </PlBlockPage>
</template>
