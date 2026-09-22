import { platforma } from "@platforma-open/milaboratories.cluster-markers.model";
import { defineAppV3 } from "@platforma-sdk/ui-vue";
import MainPage from "./pages/MainPage.vue";
// import UMAP from './pages/UMAP.vue';
import dotplot from "./pages/dotplot.vue";

export const sdkPlugin = defineAppV3(platforma, (app) => {
  return {
    progress: () => {
      return app.model.outputs.isRunning;
    },
    routes: {
      "/": () => MainPage,
      // '/umap': () => UMAP,
      "/dotplot": () => dotplot,
    },
  };
});

export const useApp = sdkPlugin.useApp;
