import { ROUTES } from "src/common/constants/routes";
import { goToPage } from "tests/utils/helpers";
import { TEST_URL } from "../common/constants";
import { getTestTag, getText } from "../utils/selectors";

describe("ToS and Privacy", () => {
  it("renders the the ToS Page", async () => {
    await goToPage(`${TEST_URL}${ROUTES.TOS}`);

    await expect(page).toHaveSelector(getText("Terms of Use"));
    await expect(page).toHaveSelector(getTestTag("czi-logo"));
    await expect(page).toHaveSelector(getTestTag("cellxgene-logo"));
  });

  it("renders the the Privacy Page", async () => {
    await goToPage(`${TEST_URL}${ROUTES.PRIVACY}`);

    await expect(page).toHaveSelector(getText("Privacy Policy"));
    await expect(page).toHaveSelector(getTestTag("czi-logo"));
    await expect(page).toHaveSelector(getTestTag("cellxgene-logo"));
  });
});
