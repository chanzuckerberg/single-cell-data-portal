/**
 * Test suite for filter-related utils.
 */

import { DEVELOPMENT_STAGE_ONTOLOGY_VIEW } from "src/components/common/Filter/common/entities";
import {
  findOntologyAncestorIds,
  findOntologyDescendantIds,
  findOntologyNodeById,
  findOntologyParentNode,
  formatNumberToScale,
  getOntologySpeciesKey,
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
  describe("Ontology", () => {
    const ONTOLOGY_ID_HUMAN_PRENATAL = "HsapDv:0000045";
    const ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN = "HsapDv:0000002";
    const ONTOLOGY_ID_HUMAN_CARNEGIE_CS1 = "HsapDv:0000003";
    const ONTOLOGY_ID_HUMAN_CLEAVAGE_CS2 = "HsapDv:0000004";
    const ONTOLOGY_ID_HUMAN_BLASTULA_CS3_5 = "HsapDv:0000006";
    const ONTOLOGY_ID_HUMAN_GASTRULA_CS6 = "HsapDv:0000010";
    const ONTOLOGY_ID_HUMAN_NEURULA_CS7_8 = "HsapDv:0000012";
    const ONTOLOGY_ID_HUMAN_ORGANOGENESIS_CS9_23 = "HsapDv:0000015";
    const ONTOLOGY_ID_HUMAN_FETAL = "HsapDv:0000037";
    const ONTOLOGY_ID_HUMAN_HUMAN_ADULT = "HsapDv:0000087";
    const ONTOLOGY_ID_HUMAN_MATURE = "HsapDv:0000204";
    const ONTOLOGY_ID_HUMAN_HUMAN_EARLY_ADULT = "HsapDv:0000088";

    const ONTOLOGY_ID_MOUSE_PRENATAL = "MmusDv:0000042";
    const ONTOLOGY_ID_OTHER_EMBRYO = "UBERON:0000068";

    describe("getOntologyKey", () => {
      it(`returns HsapDv for ${ONTOLOGY_ID_HUMAN_PRENATAL}`, () => {
        const key = "HsapDv";
        const actual = getOntologySpeciesKey(`${key}:000045`);
        expect(actual).toEqual(key);
      });
      it(`returns MmusDv for ${ONTOLOGY_ID_MOUSE_PRENATAL}`, () => {
        const key = "MmusDv";
        const actual = getOntologySpeciesKey(`${key}:000045`);
        expect(actual).toEqual(key);
      });
      it(`returns UBERON for ${ONTOLOGY_ID_OTHER_EMBRYO}`, () => {
        const key = "UBERON";
        const actual = getOntologySpeciesKey(`${key}:000045`);
        expect(actual).toEqual(key);
      });
    });
    describe("findOntologyNodeById", () => {
      [
        ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN,
        ONTOLOGY_ID_HUMAN_CARNEGIE_CS1,
        ONTOLOGY_ID_HUMAN_FETAL,
        ONTOLOGY_ID_HUMAN_HUMAN_ADULT,
      ].forEach((ontologyId) => {
        it(`finds ontology node with ID ${ontologyId}`, () => {
          const ontologyKey = getOntologySpeciesKey(ontologyId);
          const ontologyNode = findOntologyNodeById(
            DEVELOPMENT_STAGE_ONTOLOGY_VIEW[ontologyKey],
            ontologyId
          );
          expect(ontologyNode).toBeTruthy();
          // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
          expect(ontologyNode!.ontology_term_id).toEqual(ontologyId);
        });
      });
    });
    describe("findOntologyDescendants", () => {
      it(`finds descendants of ${ONTOLOGY_ID_HUMAN_PRENATAL}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_PRENATAL;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const ontologyNode = findOntologyNodeById(
          DEVELOPMENT_STAGE_ONTOLOGY_VIEW[ontologyKey],
          ontologyId
        );
        expect(ontologyNode).toBeTruthy();
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        const descendants = findOntologyDescendantIds(ontologyNode!);
        expect(descendants.length).toEqual(8);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN)
        ).toBeGreaterThanOrEqual(0);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_CARNEGIE_CS1)
        ).toBeGreaterThanOrEqual(0);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_CLEAVAGE_CS2)
        ).toBeGreaterThanOrEqual(0);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_BLASTULA_CS3_5)
        ).toBeGreaterThanOrEqual(0);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_GASTRULA_CS6)
        ).toBeGreaterThanOrEqual(0);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_NEURULA_CS7_8)
        ).toBeGreaterThanOrEqual(0);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_ORGANOGENESIS_CS9_23)
        ).toBeGreaterThanOrEqual(0);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_FETAL)
        ).toBeGreaterThanOrEqual(0);
      });
      it(`finds descendants of ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const ontologyNode = findOntologyNodeById(
          DEVELOPMENT_STAGE_ONTOLOGY_VIEW[ontologyKey],
          ontologyId
        );
        expect(ontologyNode).toBeTruthy();
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        const descendants = findOntologyDescendantIds(ontologyNode!);
        expect(descendants.length).toEqual(6);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_CARNEGIE_CS1)
        ).toBeGreaterThanOrEqual(0);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_CLEAVAGE_CS2)
        ).toBeGreaterThanOrEqual(0);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_BLASTULA_CS3_5)
        ).toBeGreaterThanOrEqual(0);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_GASTRULA_CS6)
        ).toBeGreaterThanOrEqual(0);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_NEURULA_CS7_8)
        ).toBeGreaterThanOrEqual(0);
        expect(
          descendants.indexOf(ONTOLOGY_ID_HUMAN_ORGANOGENESIS_CS9_23)
        ).toBeGreaterThanOrEqual(0);
      });
      it(`finds descendants of ${ONTOLOGY_ID_HUMAN_FETAL}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_FETAL;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const ontologyNode = findOntologyNodeById(
          DEVELOPMENT_STAGE_ONTOLOGY_VIEW[ontologyKey],
          ontologyId
        );
        expect(ontologyNode).toBeTruthy();
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        const descendants = findOntologyDescendantIds(ontologyNode!);
        expect(descendants.length).toEqual(0);
      });
    });
    describe("findOntologyAncestors", () => {
      it(`finds ancestors of ${ONTOLOGY_ID_HUMAN_CARNEGIE_CS1}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_CARNEGIE_CS1;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const rootNodes = DEVELOPMENT_STAGE_ONTOLOGY_VIEW[ontologyKey];
        const ontologyNode = findOntologyNodeById(rootNodes, ontologyId);
        expect(ontologyNode).toBeTruthy();
        const ancestorIds = [] as string[];
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        findOntologyAncestorIds(rootNodes, ontologyNode!, ancestorIds);
        expect(ancestorIds.length).toEqual(2);
        expect(
          ancestorIds.includes(ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN)
        ).toBeTruthy();
        expect(ancestorIds.includes(ONTOLOGY_ID_HUMAN_PRENATAL)).toBeTruthy();
      });
      it(`finds ancestors of ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const rootNodes = DEVELOPMENT_STAGE_ONTOLOGY_VIEW[ontologyKey];
        const ontologyNode = findOntologyNodeById(rootNodes, ontologyId);
        expect(ontologyNode).toBeTruthy();
        const ancestorIds = [] as string[];
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        findOntologyAncestorIds(rootNodes, ontologyNode!, ancestorIds);
        expect(ancestorIds.length).toEqual(1);
        expect(ancestorIds.includes(ONTOLOGY_ID_HUMAN_PRENATAL)).toBeTruthy();
      });
      it(`does not find ancestors of root node ${ONTOLOGY_ID_HUMAN_PRENATAL}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_PRENATAL;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const rootNodes = DEVELOPMENT_STAGE_ONTOLOGY_VIEW[ontologyKey];
        const ontologyNode = findOntologyNodeById(rootNodes, ontologyId);
        expect(ontologyNode).toBeTruthy();
        const ancestorIds = [] as string[];
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        findOntologyAncestorIds(rootNodes, ontologyNode!, ancestorIds);
        expect(ancestorIds.length).toEqual(0);
      });
      it(`finds ancestors of ${ONTOLOGY_ID_HUMAN_HUMAN_ADULT}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_HUMAN_ADULT;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const rootNodes = DEVELOPMENT_STAGE_ONTOLOGY_VIEW[ontologyKey];
        const ontologyNode = findOntologyNodeById(rootNodes, ontologyId);
        expect(ontologyNode).toBeTruthy();
        const ancestorIds = [] as string[];
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        findOntologyAncestorIds(rootNodes, ontologyNode!, ancestorIds);
        expect(ancestorIds.length).toEqual(1);
        expect(ancestorIds.includes(ONTOLOGY_ID_HUMAN_MATURE)).toBeTruthy();
      });
      it(`finds ancestors of ${ONTOLOGY_ID_HUMAN_HUMAN_EARLY_ADULT}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_HUMAN_EARLY_ADULT;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const rootNodes = DEVELOPMENT_STAGE_ONTOLOGY_VIEW[ontologyKey];
        const ontologyNode = findOntologyNodeById(rootNodes, ontologyId);
        expect(ontologyNode).toBeTruthy();
        const ancestorIds = [] as string[];
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        findOntologyAncestorIds(rootNodes, ontologyNode!, ancestorIds);
        expect(ancestorIds.length).toEqual(2);
        expect(
          ancestorIds.includes(ONTOLOGY_ID_HUMAN_HUMAN_ADULT)
        ).toBeTruthy();
        expect(ancestorIds.includes(ONTOLOGY_ID_HUMAN_MATURE)).toBeTruthy();
      });
    });
    describe("findOntologyParentNode", () => {
      it(`finds parent of ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const ontologyRootNodes = DEVELOPMENT_STAGE_ONTOLOGY_VIEW[ontologyKey];
        const ontologyNode = findOntologyNodeById(
          ontologyRootNodes,
          ontologyId
        );
        expect(ontologyNode).toBeTruthy();
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        const parent = findOntologyParentNode(ontologyRootNodes, ontologyNode!);
        expect(parent).toBeTruthy();
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        expect(parent!.ontology_term_id).toEqual(ONTOLOGY_ID_HUMAN_PRENATAL);
      });
      it(`finds parent of ${ONTOLOGY_ID_HUMAN_CARNEGIE_CS1}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_CARNEGIE_CS1;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const ontologyRootNodes = DEVELOPMENT_STAGE_ONTOLOGY_VIEW[ontologyKey];
        const ontologyNode = findOntologyNodeById(
          ontologyRootNodes,
          ontologyId
        );
        expect(ontologyNode).toBeTruthy();
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        const parent = findOntologyParentNode(ontologyRootNodes, ontologyNode!);
        expect(parent).toBeTruthy();
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        expect(parent!.ontology_term_id).toEqual(
          ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN
        );
      });
      it(`doesn't find parent of ${ONTOLOGY_ID_HUMAN_PRENATAL}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_PRENATAL;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const ontologyRootNodes = DEVELOPMENT_STAGE_ONTOLOGY_VIEW[ontologyKey];
        const ontologyNode = findOntologyNodeById(
          ontologyRootNodes,
          ontologyId
        );
        expect(ontologyNode).toBeTruthy();
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        const parent = findOntologyParentNode(ontologyRootNodes, ontologyNode!);
        expect(parent).toBeFalsy();
      });
    });
  });
});
