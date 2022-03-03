/**
 * Test suite for filter-related utils.
 */

import {
  DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR,
  formatNumberToScale,
  SYMBOL_MILLION,
  SYMBOL_THOUSAND,
} from "src/components/common/Filter/common/utils";

describe("filter", () => {
  describe("Format Number to Scale", () => {
    const testConfigs = [
      {
        formatted: "4",
        value: 4,
      },
      {
        formatted: "6",
        value: 6,
      },
      {
        formatted: "40",
        value: 40,
      },
      {
        formatted: "60",
        value: 60,
      },
      {
        formatted: "140",
        value: 140,
      },
      {
        formatted: "150",
        value: 150,
      },
      {
        formatted: "160",
        value: 160,
      },
      {
        formatted: "999",
        value: 999,
      },
      {
        formatted: `1.4${SYMBOL_THOUSAND}`,
        value: 1400,
      },
      {
        formatted: `1.5${SYMBOL_THOUSAND}`,
        value: 1500,
      },
      {
        formatted: `1.6${SYMBOL_THOUSAND}`,
        value: 1600,
      },
      {
        formatted: `2${SYMBOL_THOUSAND}`,
        value: 1995,
      },
      {
        formatted: `14${SYMBOL_THOUSAND}`,
        value: 14000,
      },
      {
        formatted: `15${SYMBOL_THOUSAND}`,
        value: 14500,
      },
      {
        formatted: `82${SYMBOL_THOUSAND}`,
        value: 81572,
      },
      {
        formatted: `140${SYMBOL_THOUSAND}`,
        value: 140000,
      },
      {
        formatted: `160${SYMBOL_THOUSAND}`,
        value: 160000,
      },
      {
        formatted: `190${SYMBOL_THOUSAND}`,
        value: 190000,
      },
      {
        formatted: `240${SYMBOL_THOUSAND}`,
        value: 240000,
      },
      {
        formatted: `1.4${SYMBOL_MILLION}`,
        value: 1400000,
      },
      {
        formatted: `1.6${SYMBOL_MILLION}`,
        value: 1600000,
      },
      {
        formatted: `1.9${SYMBOL_MILLION}`,
        value: 1940000,
      },
      {
        formatted: `2${SYMBOL_MILLION}`,
        value: 1950000,
      },
    ];
    testConfigs.forEach(
      ({ value, formatted }: { value: number; formatted: string }) => {
        it(`formats ${value} to ${formatted}`, () => {
          const actual = formatNumberToScale(value);
          expect(actual).toEqual(formatted);
        });
      }
    );
  });
  describe("Development Stage Leaf Ontology IDs", () => {
    it("finds leaf nodes for HsapDv:0000045", () => {
      const ancestorOntologyId = "HsapDv:0000045";
      const leafOntologyIds =
        DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ancestorOntologyId);
      expect(leafOntologyIds).toBeTruthy();
      expect(leafOntologyIds?.size).toEqual(7);
      expect(leafOntologyIds?.has("HsapDv:0000003")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000004")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000006")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000010")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000012")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000015")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000037")).toBeTruthy();
    });
    it("finds leaf nodes for HsapDv:0000002", () => {
      const ancestorOntologyId = "HsapDv:0000002";
      const leafOntologyIds =
        DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ancestorOntologyId);
      expect(leafOntologyIds).toBeTruthy();
      expect(leafOntologyIds?.size).toEqual(6);
      expect(leafOntologyIds?.has("HsapDv:0000003")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000004")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000006")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000010")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000012")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000015")).toBeTruthy();
    });
    it("finds leaf nodes for HsapDv:0000080", () => {
      const ancestorOntologyId = "HsapDv:0000080";
      const leafOntologyIds =
        DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ancestorOntologyId);
      expect(leafOntologyIds).toBeTruthy();
      expect(leafOntologyIds?.size).toEqual(3);
      expect(leafOntologyIds?.has("HsapDv:0000082")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000083")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000081")).toBeTruthy();
    });
    it("finds leaf nodes for HsapDv:0000204", () => {
      const ancestorOntologyId = "HsapDv:0000204";
      const leafOntologyIds =
        DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ancestorOntologyId);
      expect(leafOntologyIds).toBeTruthy();
      expect(leafOntologyIds?.size).toEqual(3);
      expect(leafOntologyIds?.has("HsapDv:0000086")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000088")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000091")).toBeTruthy();
    });
    it("finds leaf nodes for HsapDv:0000087", () => {
      const ancestorOntologyId = "HsapDv:0000087";
      const leafOntologyIds =
        DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ancestorOntologyId);
      expect(leafOntologyIds).toBeTruthy();
      expect(leafOntologyIds?.size).toEqual(2);
      expect(leafOntologyIds?.has("HsapDv:0000088")).toBeTruthy();
      expect(leafOntologyIds?.has("HsapDv:0000091")).toBeTruthy();
    });
    it("finds leaf nodes for MmusDv:0000042", () => {
      const ancestorOntologyId = "MmusDv:0000042";
      const leafOntologyIds =
        DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ancestorOntologyId);
      expect(leafOntologyIds).toBeTruthy();
      expect(leafOntologyIds?.size).toEqual(7);
      expect(leafOntologyIds?.has("MmusDv:0000003")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000004")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000007")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000013")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000017")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000018")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000031")).toBeTruthy();
    });
    it("finds leaf nodes for MmusDv:0000002", () => {
      const ancestorOntologyId = "MmusDv:0000002";
      const leafOntologyIds =
        DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ancestorOntologyId);
      expect(leafOntologyIds).toBeTruthy();
      expect(leafOntologyIds?.size).toEqual(6);
      expect(leafOntologyIds?.has("MmusDv:0000003")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000004")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000007")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000013")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000017")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000018")).toBeTruthy();
    });
    it("finds leaf nodes for MmusDv:0000092", () => {
      const ancestorOntologyId = "MmusDv:0000092";
      const leafOntologyIds =
        DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ancestorOntologyId);
      expect(leafOntologyIds).toBeTruthy();
      expect(leafOntologyIds?.size).toEqual(4);
      expect(leafOntologyIds?.has("MmusDv:0000036")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000112")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000061")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000097")).toBeTruthy();
    });
    it("finds leaf nodes for MmusDv:0000043", () => {
      const ancestorOntologyId = "MmusDv:0000043";
      const leafOntologyIds =
        DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ancestorOntologyId);
      expect(leafOntologyIds).toBeTruthy();
      expect(leafOntologyIds?.size).toEqual(2);
      expect(leafOntologyIds?.has("MmusDv:0000036")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000112")).toBeTruthy();
    });
    it("finds leaf nodes for MmusDv:0000110", () => {
      const ancestorOntologyId = "MmusDv:0000110";
      const leafOntologyIds =
        DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ancestorOntologyId);
      expect(leafOntologyIds).toBeTruthy();
      expect(leafOntologyIds?.size).toEqual(2);
      expect(leafOntologyIds?.has("MmusDv:0000061")).toBeTruthy();
      expect(leafOntologyIds?.has("MmusDv:0000097")).toBeTruthy();
    });
    it("finds leaf nodes for UBERON:0000068", () => {
      const ancestorOntologyId = "UBERON:0000068";
      const leafOntologyIds =
        DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ancestorOntologyId);
      expect(leafOntologyIds).toBeTruthy();
      expect(leafOntologyIds?.size).toEqual(7);
      expect(leafOntologyIds?.has("UBERON:0000106")).toBeTruthy();
      expect(leafOntologyIds?.has("UBERON:0000107")).toBeTruthy();
      expect(leafOntologyIds?.has("UBERON:0000108")).toBeTruthy();
      expect(leafOntologyIds?.has("UBERON:0000109")).toBeTruthy();
      expect(leafOntologyIds?.has("UBERON:0000110")).toBeTruthy();
      expect(leafOntologyIds?.has("UBERON:0000111")).toBeTruthy();
      expect(leafOntologyIds?.has("UBERON:0007220")).toBeTruthy();
    });
    it("finds leaf nodes for UBERON:0000092", () => {
      const ancestorOntologyId = "UBERON:0000092";
      const leafOntologyIds =
        DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ancestorOntologyId);
      expect(leafOntologyIds).toBeTruthy();
      expect(leafOntologyIds?.size).toEqual(5);
      expect(leafOntologyIds?.has("UBERON:0000069")).toBeTruthy();
      expect(leafOntologyIds?.has("UBERON:0000070")).toBeTruthy();
      expect(leafOntologyIds?.has("UBERON:0018685")).toBeTruthy();
      expect(leafOntologyIds?.has("UBERON:0000112")).toBeTruthy();
      expect(leafOntologyIds?.has("UBERON:0000113")).toBeTruthy();
    });
    it("finds leaf nodes for UBERON:0000066", () => {
      const ancestorOntologyId = "UBERON:0000066";
      const leafOntologyIds =
        DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ancestorOntologyId);
      expect(leafOntologyIds).toBeTruthy();
      expect(leafOntologyIds?.size).toEqual(2);
      expect(leafOntologyIds?.has("UBERON:0000112")).toBeTruthy();
      expect(leafOntologyIds?.has("UBERON:0000113")).toBeTruthy();
    });
  });
});
