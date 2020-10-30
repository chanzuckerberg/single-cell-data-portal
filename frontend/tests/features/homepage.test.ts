import { PROMPT_TEXT } from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/Details";
import { goToPage } from "tests/utils/helpers";
import { TEST_URL } from "../common/constants";
import { getTestTag, getText } from "../utils/selectors";

const DATASET_ROW_DOWNLOAD_BUTTON_ID = "dataset-download-button";

describe("Homepage", () => {
  it("renders the expected elements", async () => {
    await goToPage(TEST_URL);

    await expect(page).toHaveSelector(getText("cellxgene@chanzuckerberg.com"));
    await expect(page).toHaveSelector(getTestTag("logo"));
    await expect(page).toHaveSelector(getText("Dataset name"));
    await expect(page).toHaveSelector(getText("View in cellxgene"));
    await expect(page).toHaveSelector(getText("Download dataset"));
    await expect(page).toHaveSelector(getText("More information"));

    await expect(page).toHaveSelector(getTestTag("dataset-name"));
    await expect(page).toHaveSelector(getTestTag("view-dataset-link"));
    await expect(page).toHaveSelector(
      getTestTag(DATASET_ROW_DOWNLOAD_BUTTON_ID)
    );
  });

  describe("renders the download dataset modal", () => {
    it("renders the default content", async () => {
      await goToPage(TEST_URL);

      await page.click(getTestTag(DATASET_ROW_DOWNLOAD_BUTTON_ID));

      await expect(page).toHaveSelector(getText("Download Dataset"));

      await expect(page).toHaveSelector(getText("NAME"));
      expect(
        await page.innerText(getTestTag("download-asset-name"))
      ).toBeTruthy();

      await expect(page).toHaveSelector(getText("DATA FORMAT"));
      await expect(page).toHaveSelector(getText(".h5ad (AnnData v0.7)"));
      await expect(page).toHaveSelector(getText(".loom"));
      await expect(page).toHaveSelector(getText(".rds (Seurat v3)"));

      await expect(page).toHaveSelector(getText("DOWNLOAD DETAILS"));
      await expect(page).toHaveSelector(getText(PROMPT_TEXT));
    });

    it("downloads a file", async () => {
      await goToPage(TEST_URL);

      await page.click(getTestTag(DATASET_ROW_DOWNLOAD_BUTTON_ID));

      await page.click(getText(".h5ad (AnnData v0.7)"));

      const downloadLink = await page.getAttribute(
        getTestTag("download-asset-download-button"),
        "href"
      );

      expect(downloadLink).toBeTruthy();
    });
  });
});
