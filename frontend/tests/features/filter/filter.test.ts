/**
 * Test suite for filter-related queries.
 */

// App dependencies
import { expect } from "@playwright/test";
import { PublisherMetadata } from "src/common/entities";
import {
  buildSummaryCitation,
  calculateMonthsSincePublication,
  calculateRecency,
  CollectionResponse,
  createPublicationDateValues,
} from "src/common/queries/filter";
import { PUBLICATION_DATE_VALUES } from "src/components/common/Filter/common/constants";
import { test } from "tests/common/test";

const { describe } = test;

describe("filter", () => {
  describe("Calculate Months Since Publication", () => {
    test("calculates one month since publication", () => {
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
    test("calculates six months since publication", () => {
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
    test("calculates thirteen months since publication", () => {
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
  test("calculates 23 months since publication", () => {
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
    test("calculates bins for 1 month since publication", () => {
      const dateBins = createPublicationDateValues(1);
      // Expecting all date ranges (that is, 1, 3, 6, 12, 24 and 36).
      expect(dateBins.length).toEqual(6);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES);
    });
    test("calculates bins for 2 months since publication", () => {
      const dateBins = createPublicationDateValues(2);
      // Expecting 3, 6, 12, 24 and 36.
      expect(dateBins.length).toEqual(5);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES.slice(1));
    });
    test("calculates bins for 4 months since publication", () => {
      const dateBins = createPublicationDateValues(4);
      // Expecting 6, 12, 24 and 36.
      expect(dateBins.length).toEqual(4);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES.slice(2));
    });
    test("calculates bins for 7 months since publication", () => {
      const dateBins = createPublicationDateValues(7);
      // Expecting 12, 24 and 36.
      expect(dateBins.length).toEqual(3);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES.slice(3));
    });
    test("calculates bins for 32 months since publication", () => {
      const dateBins = createPublicationDateValues(32);
      // Expecting 36.
      expect(dateBins.length).toEqual(1);
      validateDateValue(dateBins, PUBLICATION_DATE_VALUES.slice(5));
    });
  });
  describe("Calculate Recency", () => {
    test("calculates recency for collection with publisher metadata", () => {
      const publishedAt = 1646092800;
      const collection = {
        publisher_metadata: {
          published_at: publishedAt,
          published_day: 11,
          published_month: 1,
          published_year: 2022,
        },
      } as CollectionResponse;
      const recency = calculateRecency(
        collection,
        collection.publisher_metadata
      );
      expect(recency).toEqual(publishedAt);
    });
    test("calculates recency for collection with revised at", () => {
      const revisedAt = 1644527777.095609; // JS month
      const collection = {
        // No publisher metadata
        revised_at: revisedAt,
      } as CollectionResponse;
      const recency = calculateRecency(
        collection,
        collection.publisher_metadata
      );
      expect(recency).toEqual(revisedAt);
    });
    test("calculates recency for collection with published at", () => {
      const publishedAt = 1644526776.095609;
      const collection = {
        // No publisher metadata or revised_at
        published_at: publishedAt,
      } as CollectionResponse;
      const recency = calculateRecency(
        collection,
        collection.publisher_metadata
      );
      expect(recency).toEqual(publishedAt);
    });
  });
  describe("Build Summary Citation", () => {
    test("builds summary citation for consortium first author", () => {
      const consortiumName = "The Tabula Sapiens Consortium";
      const journal = "bioRxiv";
      const year = 2022;
      const publisherMetadata = {
        authors: [
          { name: consortiumName },
          { family: "Quake", given: "Stephen R" },
        ],
        journal: journal,
        published_year: year,
      } as PublisherMetadata;
      const summaryCitation = buildSummaryCitation(publisherMetadata);
      expect(summaryCitation).toEqual(
        `${consortiumName} et al. (${year}) ${journal}`
      );
    });
    test("builds summary citation for person first author", () => {
      const family = "Quake";
      const journal = "bioRxiv";
      const year = 2022;
      const publisherMetadata = {
        authors: [
          { family: family, given: "Stephen R" },
          { name: "The Tabula Sapiens Consortium" },
        ],
        journal: journal,
        published_year: year,
      } as PublisherMetadata;
      const summaryCitation = buildSummaryCitation(publisherMetadata);
      expect(summaryCitation).toEqual(`${family} et al. (${year}) ${journal}`);
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
