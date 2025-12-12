import { expect, Page } from "@playwright/test";
import { TEST_URL } from "tests/common/constants";
import { goToPage, tryUntil } from "tests/utils/helpers";
import { test } from "tests/common/test";

const { describe } = test;

const GENE_EXPRESSION_SOURCE_DATA_URL = `${TEST_URL}/docs/04__Analyze%20Public%20Data/4_2__Gene%20Expression%20Documentation/4_2_6__Gene%20Expression%20Source%20Data`;
const CLOUDFRONT_BASE_URL = "https://ge-data.cellxgene.cziscience.com";
const MANIFEST_FILE = "expression-summary-files.json";

describe("Gene Expression Source Data - S3Content Component", () => {
  test("Should render S3Content component with files from manifest", async ({
    page,
  }) => {
    // Mock the manifest response
    const mockFiles = [
      "full_dataset_v1.0.0.csv",
      "condensed_dataset_v1.0.0.csv",
      "metadata_v1.0.0.json",
    ];

    await mockManifestResponse(page, mockFiles);

    await goToPage(GENE_EXPRESSION_SOURCE_DATA_URL, page);

    // Wait for the table to be visible
    await tryUntil(
      async () => {
        const table = page.locator("table");
        await expect(table).toBeVisible();
      },
      { page }
    );

    // Verify all files are listed
    for (const filename of mockFiles) {
      await expect(page.getByRole("link", { name: filename })).toBeVisible();
    }

    // Verify download links have correct href
    const firstLink = page.getByRole("link", { name: mockFiles[0] });
    await expect(firstLink).toHaveAttribute(
      "href",
      `${CLOUDFRONT_BASE_URL}/${mockFiles[0]}`
    );

    // Verify download attribute is set with filename
    await expect(firstLink).toHaveAttribute("download", mockFiles[0]);

    // Verify links open in new tab
    await expect(firstLink).toHaveAttribute("target", "_blank");
    await expect(firstLink).toHaveAttribute("rel", "noopener noreferrer");
  });

  test("Should show loading state initially", async ({ page }) => {
    // Mock a delayed response to test loading state
    await page.route(
      `${CLOUDFRONT_BASE_URL}/${MANIFEST_FILE}`,
      async (route) => {
        // Delay the response
        await new Promise((resolve) => setTimeout(resolve, 1000));
        await route.fulfill({
          status: 200,
          contentType: "application/json",
          body: JSON.stringify({ files: [] }),
        });
      }
    );

    await goToPage(GENE_EXPRESSION_SOURCE_DATA_URL, page);

    // Check that loading text appears
    await expect(page.getByText("Loading file listing...")).toBeVisible();
  });

  test("Should show error message when manifest fetch fails", async ({
    page,
  }) => {
    // Mock a failed response
    await page.route(
      `${CLOUDFRONT_BASE_URL}/${MANIFEST_FILE}`,
      async (route) => {
        await route.fulfill({
          status: 404,
          contentType: "text/plain",
          body: "Not Found",
        });
      }
    );

    await goToPage(GENE_EXPRESSION_SOURCE_DATA_URL, page);

    // Wait for error message to appear
    await tryUntil(
      async () => {
        await expect(
          page.getByText(/Error:.*Failed to fetch manifest/)
        ).toBeVisible();
      },
      { page }
    );
  });

  test("Should show error message when manifest is not valid JSON", async ({
    page,
  }) => {
    // Mock an invalid JSON response
    await page.route(
      `${CLOUDFRONT_BASE_URL}/${MANIFEST_FILE}`,
      async (route) => {
        await route.fulfill({
          status: 200,
          contentType: "application/json",
          body: "invalid json {",
        });
      }
    );

    await goToPage(GENE_EXPRESSION_SOURCE_DATA_URL, page);

    // Wait for error message to appear
    await tryUntil(
      async () => {
        await expect(
          page.getByText(/Error:.*Manifest file is not valid JSON/)
        ).toBeVisible();
      },
      { page }
    );
  });

  test("Should show 'No files found' when files array is empty", async ({
    page,
  }) => {
    await mockManifestResponse(page, []);

    await goToPage(GENE_EXPRESSION_SOURCE_DATA_URL, page);

    // Wait for no files message to appear
    await tryUntil(
      async () => {
        await expect(page.getByText("No files found")).toBeVisible();
      },
      { page }
    );
  });

  test("Should handle manifest without files property", async ({ page }) => {
    // Mock a response without files property
    await page.route(
      `${CLOUDFRONT_BASE_URL}/${MANIFEST_FILE}`,
      async (route) => {
        await route.fulfill({
          status: 200,
          contentType: "application/json",
          body: JSON.stringify({ other: "data" }),
        });
      }
    );

    await goToPage(GENE_EXPRESSION_SOURCE_DATA_URL, page);

    // Should show no files found
    await tryUntil(
      async () => {
        await expect(page.getByText("No files found")).toBeVisible();
      },
      { page }
    );
  });

  test("Should render documentation page content", async ({ page }) => {
    await mockManifestResponse(page, ["test-file.csv"]);

    await goToPage(GENE_EXPRESSION_SOURCE_DATA_URL, page);

    // Verify page heading
    await expect(
      page.getByRole("heading", { name: "Gene Expression Source Data" })
    ).toBeVisible();

    // Verify description text is present
    await expect(page.getByText(/Full dataset/)).toBeVisible();
    await expect(page.getByText(/Condensed dataset/)).toBeVisible();
  });
});

/**
 * Helper function to mock the manifest response
 */
async function mockManifestResponse(page: Page, files: string[]) {
  await page.route(`${CLOUDFRONT_BASE_URL}/${MANIFEST_FILE}`, async (route) => {
    await route.fulfill({
      status: 200,
      contentType: "application/json",
      body: JSON.stringify({ files }),
    });
  });
}
