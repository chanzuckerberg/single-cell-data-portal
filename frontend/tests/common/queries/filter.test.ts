/**
 * Test suite for filter-related queries.
 */

// App dependencies
import {
  calculateMonthsSincePublication,
  createPublicationDateValues,
} from "src/common/queries/filter";
import { PUBLICATION_DATE_VALUES } from "src/components/common/Filter/common/entities";

describe("filter", () => {
  describe("Calculate Months Since Publication", () => {
    it("calculates one month since publication", () => {
      // Publication date - 12/2021
      // Today's date - 1/2022
      const monthsSincePublication = calculateMonthsSincePublication(
        1,
        2022,
        12,
        2021
      );
      expect(monthsSincePublication).toEqual(1);
    });
    it("calculates six months since publication", () => {
      // Publication date - 12/2021
      // Today's date - 1/2022
      const monthsSincePublication = calculateMonthsSincePublication(
        1,
        2022,
        7,
        2021
      );
      expect(monthsSincePublication).toEqual(6);
    });
    it("calculates thirteen months since publication", () => {
      // Publication date - 12/2021
      // Today's date - 1/2022
      const monthsSincePublication = calculateMonthsSincePublication(
        1,
        2022,
        12,
        2020
      );
      expect(monthsSincePublication).toEqual(13);
    });
  });
  it("calculates 23 months since publication", () => {
    // Publication date - 12/2021
    // Today's date - 1/2022
    const monthsSincePublication = calculateMonthsSincePublication(
      12,
      2022,
      1,
      2021
    );
    expect(monthsSincePublication).toEqual(23);
  });
  describe("Calculate Date Bins", () => {
    it("calculates bins for 1 month since publication", () => {
      const dateBins = createPublicationDateValues(1);
      // Expecting all date ranges (that is, 1, 3, 6, 12, 24 and 36).
      expect(dateBins.length).toEqual(6);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES);
    });
    it("calculates bins for 2 months since publication", () => {
      const dateBins = createPublicationDateValues(2);
      // Expecting 3, 6, 12, 24 and 36.
      expect(dateBins.length).toEqual(5);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES.slice(1));
    });
    it("calculates bins for 4 months since publication", () => {
      const dateBins = createPublicationDateValues(4);
      // Expecting 6, 12, 24 and 36.
      expect(dateBins.length).toEqual(4);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES.slice(2));
    });
    it("calculates bins for 7 months since publication", () => {
      const dateBins = createPublicationDateValues(7);
      // Expecting 12, 24 and 36.
      expect(dateBins.length).toEqual(3);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES.slice(3));
    });
    it("calculates bins for 32 months since publication", () => {
      const dateBins = createPublicationDateValues(32);
      // Expecting 36.
      expect(dateBins.length).toEqual(1);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES.slice(5));
    });
  });
});

/**
 * Check each returned, actual date value and confirm it is in the expected.
 * @param actualDataValues - Array of numbers (in the result) to check.
 * @param expectedDateValues - Array of numbers that are expected in the result.
 */
function validateDateValue(
  actualDataValues: number[],
  expectedDateValues: number[]
) {
  actualDataValues.forEach((dateValue: number, index: number) => {
    expect(dateValue).toEqual(expectedDateValues[index]);
  });
}
