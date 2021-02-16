import { PROMPT_TEXT } from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/Details";
import { goToPage, tryUntil } from "tests/utils/helpers";
import { getTestTag, getText } from "../utils/selectors";

const DATASET_ROW_DOWNLOAD_BUTTON_ID = "dataset-download-button";

describe("Homepage", () => {
  it("renders the expected elements", async () => {
    await goToPage();

    await expect(page).toHaveText("cellxgene@chanzuckerberg.com");
    await expect(page).toHaveSelector(getTestTag("logo"));
    await expect(page).toHaveText("Dataset name");
    await expect(page).toHaveText("View in cellxgene");
    await expect(page).toHaveText("Download dataset");
    await expect(page).toHaveText("More information");

    await expect(page).toHaveSelector(getTestTag("dataset-name"));
    await expect(page).toHaveSelector(getTestTag("view-dataset-link"));
    await expect(page).toHaveSelector(
      getTestTag(DATASET_ROW_DOWNLOAD_BUTTON_ID)
    );
  });

  describe("renders the download dataset modal", () => {
    it("renders the default content", async () => {
      await goToPage();

      await page.click(getTestTag(DATASET_ROW_DOWNLOAD_BUTTON_ID));

      await expect(page).toHaveText("Download Dataset");

      await expect(page).toHaveText("NAME");
      expect(
        await page.innerText(getTestTag("download-asset-name"))
      ).toBeTruthy();

      await expect(page).toHaveText("DATA FORMAT");
      await expect(page).toHaveText(".h5ad (AnnData v0.7)");
      await expect(page).toHaveText(".loom");
      await expect(page).toHaveText(".rds (Seurat v3)");

      await expect(page).toHaveText("DOWNLOAD DETAILS");
      await expect(page).toHaveText(PROMPT_TEXT);
    });

    it("downloads a file", async () => {
      await goToPage();

      await page.click(getTestTag(DATASET_ROW_DOWNLOAD_BUTTON_ID));

      await page.click(getText(".h5ad (AnnData v0.7)"));

      await tryUntil(async () => {
        const downloadLink = await page.getAttribute(
          getTestTag("download-asset-download-button"),
          "href"
        );

        expect(downloadLink).toBeTruthy();
      });
    });
  });
});
