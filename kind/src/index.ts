import { assertParamsObject, defineBlockKind } from "@platforma-sdk/block-kind";
import type { PlRef } from "@platforma-sdk/model";
import { isPlRef } from "@platforma-sdk/model";
import { name, version } from "../package.json" with { type: "json" };

/**
 * This block's init-params contract — the five settings a user picks by hand:
 * the upstream cluster annotation markers are called against, how many top
 * markers to keep per cluster, the two significance thresholds, and the
 * specificity mode.
 *
 * `title` is absent. It is a field of `BlockArgs` that no part of the UI ever
 * writes, so it is always undefined and there is nothing for a template to
 * carry. The per-instance label belongs in the defaultBlockLabel /
 * customBlockLabel pattern, not here.
 *
 * Every field is optional. A block with no cluster annotation chosen is
 * ordinary state the UI reaches — it is exactly what a freshly created block
 * looks like — and the projection hands that state back untouched, so a
 * required field would make this block export a file its own kind refuses.
 */
export type BlockParams = {
  clusterAnnotationRef?: PlRef;
  topN?: number;
  logfcCutoff?: number;
  pvalCutoff?: number;
  strictOverlap?: boolean;
};

// Identity (`name`/`version`) comes from this package's own `package.json`, so
// the on-wire `{name}@{version}` reference can never drift from what npm
// publishes; the bundler inlines the JSON import.
export const kind = defineBlockKind<BlockParams>({
  name,
  version,
  parseInitializationParams,
});

// Internals

/** The same contract at runtime, for params arriving from a template file rather than typed code. */
function parseInitializationParams(value: unknown): BlockParams {
  assertParamsObject(value);

  const { clusterAnnotationRef, topN, logfcCutoff, pvalCutoff, strictOverlap } = value;

  if (clusterAnnotationRef !== undefined && !isPlRef(clusterAnnotationRef)) {
    throw new Error(
      "'clusterAnnotationRef' must be a reference to an upstream column, written as { __isRef: true, blockId, name }.",
    );
  }
  // Markers kept per cluster — a count.
  if (topN !== undefined && !isPositiveInteger(topN)) {
    throw new Error("'topN' must be an integer greater than or equal to 1.");
  }
  // A floor on |log2 fold change|, so any non-negative number.
  if (logfcCutoff !== undefined && !isNonNegativeNumber(logfcCutoff)) {
    throw new Error("'logfcCutoff' must be a number greater than or equal to 0.");
  }
  // A ceiling on the adjusted p-value, so a probability.
  if (pvalCutoff !== undefined && !isProbability(pvalCutoff)) {
    throw new Error("'pvalCutoff' must be a number between 0 and 1.");
  }
  if (strictOverlap !== undefined && typeof strictOverlap !== "boolean") {
    throw new Error("'strictOverlap' must be true or false.");
  }

  return { clusterAnnotationRef, topN, logfcCutoff, pvalCutoff, strictOverlap };
}

function isPositiveInteger(value: unknown): value is number {
  return typeof value === "number" && Number.isInteger(value) && value >= 1;
}

function isNonNegativeNumber(value: unknown): value is number {
  return typeof value === "number" && Number.isFinite(value) && value >= 0;
}

function isProbability(value: unknown): value is number {
  return typeof value === "number" && Number.isFinite(value) && value >= 0 && value <= 1;
}
