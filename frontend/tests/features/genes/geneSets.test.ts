import { expect, Page, test } from "@playwright/test";
import { goToPage } from "tests/utils/helpers";

const { describe, skip } = test;

// Excluding gene sets test suite on go live of filter (#2391 and 1718) as the original Portal homepage that contained
// the gene sets upload functionality will no longer available to users. Remove this test suite with #2121.
skip("Gene sets", () => {
  describe("CSV Validation", () => {
    describe("Given a non-UTF CSV", () => {
      test("returns an error message", async ({ page }) => {
        await assertCsvResult(
          "/fixtures/geneSets/Tidy CSV With Comments non-UTF.csv",
          page
        );
      });
    });

    describe("Given a UTF-16 CSV", () => {
      test("returns no error message and no comments in the result", async ({
        page,
      }) => {
        await assertCsvResult(
          "/fixtures/geneSets/Tidy CSV With Comments UTF16.csv",
          page
        );
      });
    });

    describe("Given a UTF-8 CSV", () => {
      test("returns no error message", async ({ page }) => {
        await assertCsvResult("/fixtures/geneSets/Tidy CSV UTF8.csv", page);
      });

      describe("with comments", () => {
        test("returns no error message and no comments in the result", async ({
          page,
        }) => {
          await assertCsvResult(
            "/fixtures/geneSets/Tidy CSV With Comments UTF8.csv",
            page
          );
        });
      });
    });

    describe("Given a UTF-8 with comments, empty gene set name, empty gene set description, empty gene symbol, empty gene description, duplicate gene set name, and duplicate gene symbol", () => {
      test("returns multiple errors", async ({ page }) => {
        await assertCsvResult(
          "/fixtures/geneSets/Tidy CSV With Comments, empty gene set name, empty gene set description, empty gene symbol, empty gene description, duplicate gene set name, and duplicate gene symbol UTF8.csv",
          page
        );
      });
    });

    describe("Given a UTF-8 with comments and existing gene set name in database", () => {
      test("returns one error", async ({ page }) => {
        await assertCsvResult(
          "/fixtures/geneSets/Tidy CSV With Comments and existing gene set name in database UTF8.csv",
          page
        );
      });
    });

    describe("Given a UTF-8 with illegal whitespace", () => {
      test("returns multiple error", async ({ page }) => {
        await assertCsvResult(
          "/fixtures/geneSets/Tidy CSV With Whitespace UTF8.csv",
          page
        );
      });
    });
  });
});

async function assertCsvResult(filename: string, page: Page) {
  await goToPage(undefined, page);

  await Promise.all([
    page.on("filechooser", async (fileChooser) => {
      fileChooser.setFiles(__dirname + filename);
    }),

    page.getByTestId("upload-csv").click(),
  ]);

  const csvResult = await page.getByTestId("csv-result").textContent();

  await expect(csvResult).toMatchSnapshot();
}
