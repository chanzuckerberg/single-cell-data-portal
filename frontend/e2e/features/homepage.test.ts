import { goToPage } from "e2e/utils/helpers";
import { TEST_URL } from "../common/constants";
import { getTestTag, getText } from "../utils/selectors";

describe("Homepage", () => {
  it("renders the expected elements", async () => {
    await goToPage(TEST_URL);

    await expect(page).toHaveSelector(getTestTag("logo"));
    await expect(page).toHaveSelector(getText("Name of dataset"));
    await expect(page).toHaveSelector(getText("View dataset in cellxgene"));
    await expect(page).toHaveSelector(getText("More information"));

    await expect(page).toHaveSelector(getTestTag("dataset-name"));
    await expect(page).toHaveSelector(getTestTag("view-dataset-link"));
  });
});
