<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAgDataTable, PlAgDataTableToolsPanel, PlBlockPage, PlBtnGhost, PlDropdownRef, PlMaskIcon24, PlSlideModal } from '@platforma-sdk/ui-vue';
import type { PlDataTableSettings } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';
import { computed, ref } from 'vue';

const app = useApp();

const tableSettings = computed<PlDataTableSettings>(() => ({
  sourceType: 'ptable',
  pTable: app.model.outputs.clusterMarkersPt?.table,
  sheets: app.model.outputs.clusterMarkersPt?.sheets,
}));

// const settingsAreShown = ref(app.model.outputs.UMAPPf === undefined)
const settingsAreShown = ref(true);
const showSettings = () => {
  settingsAreShown.value = true;
};

// const covariateOptions = computed(() => {
//   return app.model.outputs.metadataOptions?.map((v) => ({
//     value: v.ref,
//     label: v.label,
//   })) ?? [];
// });

</script>

<template>
  <PlBlockPage>
    <template #title>Cluster Markers</template>
    <template #append>
      <!-- PlAgDataTableToolsPanel controls showing  Export column and filter-->
      <PlAgDataTableToolsPanel/>
      <PlBtnGhost @click.stop="showSettings">
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
    <PlSlideModal v-model="settingsAreShown">
      <template #title>Settings</template>
      <PlDropdownRef
        v-model="app.model.args.countsRef" :options="app.model.outputs.countsOptions"
        label="Select dataset"
      />
      <PlDropdownRef
        v-model="app.model.args.clusterAnnotationRef" :options="app.model.outputs.clusterAnnotationOptions"
        label="Cluster annotation"
      />
    </PlSlideModal>
  </PlBlockPage>
</template>
