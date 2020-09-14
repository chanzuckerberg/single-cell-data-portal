import { goToPage } from "tests/utils/helpers";
import { TEST_URL } from "../common/constants";
import { getTestTag, getText } from "../utils/selectors";

describe("Homepage", () => {
  it("renders the expected elements", async () => {
    await goToPage(TEST_URL);

    await expect(page).toHaveSelector(getTestTag("logo"));
    await expect(page).toHaveSelector(getText("Dataset name"));
    await expect(page).toHaveSelector(getText("View in cellxgene"));
    await expect(page).toHaveSelector(getText("Download dataset"));
    await expect(page).toHaveSelector(getText("More information"));

    await expect(page).toHaveSelector(getTestTag("dataset-name"));
    await expect(page).toHaveSelector(getTestTag("view-dataset-link"));
    await expect(page).toHaveSelector(getTestTag("dataset-download-button"));
  });
});
