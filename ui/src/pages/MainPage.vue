<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAgDataTableV2, PlAlert, PlBlockPage, PlBtnGhost, PlBtnGroup, PlDropdownRef, PlMaskIcon24, PlNumberField, PlRow, PlSlideModal, usePlDataTableSettingsV2 } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';
import { reactive } from 'vue';
import type { PlRef } from '@platforma-sdk/model';
import { plRefsEqual } from '@platforma-sdk/model';

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
  settingsOpen: app.model.args.countsRef === undefined,
});

function setInput(inputRef?: PlRef) {
  app.model.args.countsRef = inputRef;
  if (inputRef)
    app.model.args.title = app.model.outputs.countsOptions?.find((o) => plRefsEqual(o.ref, inputRef))?.label;
  else
    app.model.args.title = undefined;
}

</script>

<template>
  <PlBlockPage>
    <template #title>Cluster Markers</template>
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
        v-model="app.model.args.countsRef" :options="app.model.outputs.countsOptions"
        label="Select dataset"
        clearable @update:model-value="setInput"
      />
      <PlDropdownRef
        v-model="app.model.args.clusterAnnotationRef" :options="app.model.outputs.clusterAnnotationOptions"
        label="Cluster annotation"
      />
      <PlNumberField
        v-model="app.model.args.topN"
        label="Number of top markers" :minValue="1" :step="1"
      >
        <template #tooltip>
          Select number of top markers to visualize.
        </template>
      </PlNumberField>
      <PlBtnGroup
        v-model="app.model.args.strictOverlap"
        label="Overlap filtering"
        :options="overlapOptions"
      >
        <template #tooltip>
          <div>
            <strong>Overlap Filtering</strong><br/>
            Controls how strictly cluster markers are filtered based on expression overlap between clusters.<br/><br/>
            <strong>Non-exclusive:</strong> Genes are considered markers if expressed in at least 20% of cells in the cluster.<br/><br/>
            <strong>Strict overlap:</strong> Genes are considered markers only if expressed in at least 20% of cells in the cluster AND less than 20% of cells in other clusters.<br/><br/>
          </div>
        </template>
      </PlBtnGroup>
      <PlRow>
        <PlNumberField
          v-model="app.model.args.logfcCutoff"
          label="Log2(FC)" :minValue="0" :step="0.1"
        >
          <template #tooltip>
            Select a valid absolute log2(FC) threshold for identifying
            significant cluster markers.
          </template>
        </PlNumberField>
        <PlNumberField
          v-model="app.model.args.pvalCutoff"
          label="Adjusted p-value" :minValue="0" :maxValue="1" :step="0.01"
        />
      </PlRow>
      <!-- Add warnings if selected threshold are out of most commonly used bounds -->
      <PlAlert v-if="app.model.args.pvalCutoff > 0.05" type="warn">
        {{ "Warning: The selected adjusted p-value threshold is higher than the most commonly used 0.05" }}
      </PlAlert>
      <PlAlert v-if="app.model.args.logfcCutoff < 0.6" type="warn">
        {{ "Warning: The selected Log2(FC) threshold may be too low for most use cases" }}
      </PlAlert>
    </PlSlideModal>
  </PlBlockPage>
</template>
