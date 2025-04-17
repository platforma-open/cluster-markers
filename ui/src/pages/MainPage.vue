<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAgDataTable, PlAgDataTableToolsPanel, PlNumberField, PlBlockPage, PlBtnGhost, PlDropdownRef, PlMaskIcon24, PlSlideModal, PlRow, PlAlert } from '@platforma-sdk/ui-vue';
import type { PlDataTableSettings } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';
import { computed, reactive } from 'vue';
import type { PlRef } from '@platforma-sdk/model';
import { plRefsEqual } from '@platforma-sdk/model';

const app = useApp();

const tableSettings = computed<PlDataTableSettings>(() => ({
  sourceType: 'ptable',
  pTable: app.model.outputs.clusterMarkersPt?.table,
  sheets: app.model.outputs.clusterMarkersPt?.sheets,
}));

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
      <!-- PlAgDataTableToolsPanel controls showing  Export column and filter-->
      <PlAgDataTableToolsPanel/>
      <PlBtnGhost @click.stop="() => data.settingsOpen = true">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>
    <PlAgDataTable
      v-model="app.model.ui.tableState"
      :settings="tableSettings"
      show-columns-panel
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
