/**
 * Test suite for filter-related utils.
 */

import {
  formatNumberToScale,
  SYMBOL_MILLION,
  SYMBOL_THOUSAND,
} from "src/components/common/Filter/common/utils";

describe("filter", () => {
  describe("Format Number to Magnitude", () => {
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
        formatted: `0.1${SYMBOL_MILLION}`,
        value: 140000,
      },
      {
        formatted: `0.2${SYMBOL_MILLION}`,
        value: 160000,
      },
      {
        formatted: `0.2${SYMBOL_MILLION}`,
        value: 190000,
      },
      {
        formatted: `0.2${SYMBOL_MILLION}`,
        value: 240000,
      },
      {
        formatted: `0.3${SYMBOL_MILLION}`,
        value: 250000,
      },
      {
        formatted: `0.3${SYMBOL_MILLION}`,
        value: 260000,
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
});
