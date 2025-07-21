import { model } from '@platforma-open/milaboratories.cluster-markers.model';
import { defineApp } from '@platforma-sdk/ui-vue';
import MainPage from './pages/MainPage.vue';
// import UMAP from './pages/UMAP.vue';
import dotplot from './pages/dotplot.vue';

export const sdkPlugin = defineApp(model, () => {
  return {
    routes: {
      '/': () => MainPage,
      // '/umap': () => UMAP,
      '/dotplot': () => dotplot,
    },
  };
});

export const useApp = sdkPlugin.useApp;
