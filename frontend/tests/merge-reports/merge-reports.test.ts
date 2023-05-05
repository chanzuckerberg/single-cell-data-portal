import { test } from "@playwright/test";
import { mergeHTMLReports } from "playwright-merge-html-reports";

mergeHTMLReports([
  process.cwd() + "/html-reports-1",
  process.cwd() + "/html-reports-2",
  process.cwd() + "/html-reports-3",
]);

const inputReportPaths = [
  process.cwd() + "/html-reports-1",
  process.cwd() + "/html-reports-2",
  process.cwd() + "/html-reports-3",
];

const config = {
  outputFolderName: "html-report", // default value
  outputBasePath: process.cwd(), // default value
};
test.describe("Merge reports", () => {
  test("Should merge Playwright shard reports", async () => {
    mergeHTMLReports(inputReportPaths, config);
  });
});
