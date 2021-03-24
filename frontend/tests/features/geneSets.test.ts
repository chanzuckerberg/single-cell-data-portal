import { goToPage } from "tests/utils/helpers";
import { getTestID } from "tests/utils/selectors";

describe("Gene sets", () => {
  describe("CSV Validation", () => {
    describe("Given a non-UTF CSV", () => {
      it("returns an error message", async () => {
        await assertCsvResult(
          "/fixtures/geneSets/Tidy CSV With Comments non-UTF.csv"
        );
      });
    });

    describe("Given a UTF-16 CSV", () => {
      it("returns no error message and no comments in the result", async () => {
        await assertCsvResult(
          "/fixtures/geneSets/Tidy CSV With Comments UTF16.csv"
        );
      });
    });

    describe("Given a UTF-8 CSV", () => {
      it("returns no error message", async () => {
        await assertCsvResult("/fixtures/geneSets/Tidy CSV UTF8.csv");
      });

      describe("with comments", () => {
        it("returns no error message and no comments in the result", async () => {
          await assertCsvResult(
            "/fixtures/geneSets/Tidy CSV With Comments UTF8.csv"
          );
        });
      });
    });

    describe("Given a UTF-8 with comments, empty gene set name, empty gene set description, empty gene symbol, empty gene description, duplicate gene set name, and duplicate gene symbol", () => {
      it("returns multiple errors", async () => {
        await assertCsvResult(
          "/fixtures/geneSets/Tidy CSV With Comments, empty gene set name, empty gene set description, empty gene symbol, empty gene description, duplicate gene set name, and duplicate gene symbol UTF8.csv"
        );
      });
    });

    describe("Given a UTF-8 with comments and existing gene set name in database", () => {
      it("returns one error", async () => {
        await assertCsvResult(
          "/fixtures/geneSets/Tidy CSV With Comments and existing gene set name in database UTF8.csv"
        );
      });
    });

    describe("Given a UTF-8 with illegal whitespace", () => {
      it("returns multiple error", async () => {
        await assertCsvResult(
          "/fixtures/geneSets/Tidy CSV With Whitespace UTF8.csv"
        );
      });
    });
  });
});

async function assertCsvResult(filename: string) {
  await goToPage();

  await Promise.all([
    page.on("filechooser", async (fileChooser) => {
      fileChooser.setFiles(__dirname + filename);
    }),

    page.click(getTestID("upload-csv")),
  ]);

  const csvResult = await page.innerText(getTestID("csv-result"));

  expect(csvResult).toMatchSnapshot();
}
