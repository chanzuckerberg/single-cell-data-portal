/**
 * Test suite for filter-related utils.
 */

import test, { expect } from "@playwright/test";
import { DEVELOPMENT_STAGE_ONTOLOGY_TERM_SET } from "src/components/common/Filter/common/constants";
import {
  findOntologyNodeById,
  findOntologyParentNode,
  formatNumberToScale,
  getOntologySpeciesKey,
  SYMBOL_MILLION,
  SYMBOL_THOUSAND,
} from "src/components/common/Filter/common/utils";

const { describe } = test;

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
        test(`formats ${value} to ${formatted}`, () => {
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
    const ONTOLOGY_ID_HUMAN_FETAL = "HsapDv:0000037";
    const ONTOLOGY_ID_HUMAN_HUMAN_ADULT = "HsapDv:0000087";

    const ONTOLOGY_ID_MOUSE_PRENATAL = "MmusDv:0000042";
    const ONTOLOGY_ID_OTHER_EMBRYO = "UBERON:0000068";

    describe("getOntologyKey", () => {
      test(`returns HsapDv for ${ONTOLOGY_ID_HUMAN_PRENATAL}`, () => {
        const key = "HsapDv";
        const actual = getOntologySpeciesKey(`${key}:000045`);
        expect(actual).toEqual(key);
      });
      test(`returns MmusDv for ${ONTOLOGY_ID_MOUSE_PRENATAL}`, () => {
        const key = "MmusDv";
        const actual = getOntologySpeciesKey(`${key}:000045`);
        expect(actual).toEqual(key);
      });
      test(`returns UBERON for ${ONTOLOGY_ID_OTHER_EMBRYO}`, () => {
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
        test(`finds ontology node with ID ${ontologyId}`, () => {
          const ontologyKey = getOntologySpeciesKey(ontologyId);
          const ontologyRootNodes =
            DEVELOPMENT_STAGE_ONTOLOGY_TERM_SET[ontologyKey];
          expect(ontologyRootNodes).toBeTruthy();
          const ontologyNode = findOntologyNodeById(
            // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
            ontologyRootNodes!,
            ontologyId
          );
          expect(ontologyNode).toBeTruthy();
          // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
          expect(ontologyNode!.ontology_term_id).toEqual(ontologyId);
        });
      });
    });
    describe("findOntologyParentNode", () => {
      test(`finds parent of ${ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const ontologyRootNodes =
          DEVELOPMENT_STAGE_ONTOLOGY_TERM_SET[ontologyKey];
        expect(ontologyRootNodes).toBeTruthy();
        const ontologyNode = findOntologyNodeById(
          // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
          ontologyRootNodes!,
          ontologyId
        );
        expect(ontologyNode).toBeTruthy();
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy checks above
        const parent = findOntologyParentNode(
          ontologyRootNodes!,
          ontologyNode!
        );
        expect(parent).toBeTruthy();
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        expect(parent!.ontology_term_id).toEqual(ONTOLOGY_ID_HUMAN_PRENATAL);
      });
      test(`finds parent of ${ONTOLOGY_ID_HUMAN_CARNEGIE_CS1}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_CARNEGIE_CS1;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const ontologyRootNodes =
          DEVELOPMENT_STAGE_ONTOLOGY_TERM_SET[ontologyKey];
        expect(ontologyRootNodes).toBeTruthy();
        const ontologyNode = findOntologyNodeById(
          // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
          ontologyRootNodes!,
          ontologyId
        );
        expect(ontologyNode).toBeTruthy();
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy checks above
        const parent = findOntologyParentNode(
          ontologyRootNodes!,
          ontologyNode!
        );
        expect(parent).toBeTruthy();
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
        expect(parent!.ontology_term_id).toEqual(
          ONTOLOGY_ID_HUMAN_EMBRYONIC_HUMAN
        );
      });
      test(`doesn't find parent of ${ONTOLOGY_ID_HUMAN_PRENATAL}`, () => {
        const ontologyId = ONTOLOGY_ID_HUMAN_PRENATAL;
        const ontologyKey = getOntologySpeciesKey(ontologyId);
        const ontologyRootNodes =
          DEVELOPMENT_STAGE_ONTOLOGY_TERM_SET[ontologyKey];
        expect(ontologyRootNodes).toBeTruthy();
        const ontologyNode = findOntologyNodeById(
          // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy check above
          ontologyRootNodes!,
          ontologyId
        );
        expect(ontologyNode).toBeTruthy();
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- truthy checks above
        const parent = findOntologyParentNode(
          ontologyRootNodes!,
          ontologyNode!
        );
        expect(parent).toBeFalsy();
      });
    });
  });
});
