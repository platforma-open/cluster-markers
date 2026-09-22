import { describe, expect, it } from "vitest";
import { kind } from "./index";

const parse = (params: unknown) => kind.parseInitializationParams(params);

/** A cluster annotation ref exactly as `resultPool.getOptions` mints it for the dropdown. */
const CLUSTER_REF = {
  __isRef: true,
  blockId: "leiden-clustering-1",
  name: "pl7.app/rna-seq/leidencluster",
} as const;

describe("the params envelope", () => {
  it.each([
    ["null", null],
    ["a number", 5],
    ["an array", ["topN"]],
    ["a string", "topN=3"],
  ])("rejects %s", (_label, value) => {
    expect(() => parse(value)).toThrow();
  });

  it("accepts an empty object — a block created with nothing seeded", () => {
    expect(parse({})).toEqual({});
  });

  it("drops keys the contract does not name", () => {
    expect(parse({ topN: 5, title: "Cluster Markers - CD8" })).toEqual({ topN: 5 });
  });
});

describe("clusterAnnotationRef", () => {
  it("accepts a PlRef", () => {
    expect(parse({ clusterAnnotationRef: CLUSTER_REF })).toEqual({
      clusterAnnotationRef: CLUSTER_REF,
    });
  });

  it.each([
    ["an object missing the ref marker", { blockId: "b1", name: "pl7.app/rna-seq/cellType" }],
    ["a bare column name", "pl7.app/rna-seq/cellType"],
    ["null", null],
  ])("rejects %s", (_label, ref) => {
    expect(() => parse({ clusterAnnotationRef: ref })).toThrow("'clusterAnnotationRef' must be");
  });
});

describe("topN", () => {
  it("accepts a positive integer", () => {
    expect(parse({ topN: 3 })).toEqual({ topN: 3 });
  });

  it.each([
    ["zero", 0],
    ["a negative count", -1],
    ["a fraction", 2.5],
    ["a numeric string", "3"],
  ])("rejects %s", (_label, topN) => {
    expect(() => parse({ topN })).toThrow("'topN' must be an integer greater than or equal to 1.");
  });
});

describe("logfcCutoff", () => {
  it.each([
    ["the default", 1],
    ["zero — no fold-change filtering", 0],
    ["a fraction", 0.6],
  ])("accepts %s", (_label, logfcCutoff) => {
    expect(parse({ logfcCutoff })).toEqual({ logfcCutoff });
  });

  it.each([
    ["a negative threshold", -1],
    ["NaN", Number.NaN],
    ["a numeric string", "1.0"],
  ])("rejects %s", (_label, logfcCutoff) => {
    expect(() => parse({ logfcCutoff })).toThrow("'logfcCutoff' must be");
  });
});

describe("pvalCutoff", () => {
  it.each([
    ["the default", 0.01],
    ["the boundary", 1],
  ])("accepts %s", (_label, pvalCutoff) => {
    expect(parse({ pvalCutoff })).toEqual({ pvalCutoff });
  });

  it.each([
    ["a value above 1", 1.5],
    ["a negative value", -0.01],
    ["a numeric string", "0.05"],
  ])("rejects %s", (_label, pvalCutoff) => {
    expect(() => parse({ pvalCutoff })).toThrow("'pvalCutoff' must be a number between 0 and 1.");
  });
});

describe("strictOverlap", () => {
  it.each([
    ["non-exclusive", false],
    ["strict overlap", true],
  ])("accepts %s", (_label, strictOverlap) => {
    expect(parse({ strictOverlap })).toEqual({ strictOverlap });
  });

  it("rejects a string", () => {
    expect(() => parse({ strictOverlap: "true" })).toThrow(
      "'strictOverlap' must be true or false.",
    );
  });
});

it("round-trips a fully configured block", () => {
  const params = {
    clusterAnnotationRef: CLUSTER_REF,
    topN: 10,
    logfcCutoff: 0.8,
    pvalCutoff: 0.05,
    strictOverlap: true,
  };

  expect(parse(params)).toEqual(params);
});
