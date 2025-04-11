<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAgDataTable, PlAgDataTableToolsPanel, PlBlockPage, PlBtnGhost, PlDropdownRef, PlMaskIcon24, PlSlideModal } from '@platforma-sdk/ui-vue';
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
    </PlSlideModal>
  </PlBlockPage>
</template>
