import { ROUTES } from "src/common/constants/routes";
import { goToPage } from "tests/utils/helpers";
import { TEST_URL } from "../common/constants";
import { getTestID, getText } from "../utils/selectors";

describe("ToS and Privacy", () => {
  it("renders the the ToS Page", async () => {
    await goToPage(`${TEST_URL}${ROUTES.TOS}`);
    page.screenshot({ path: `./tmp-screenshots/tos-${Date.now()}.png` });

    await expect(page).toHaveSelector(getText("Terms of Use"));
    await expect(page).toHaveSelector(getTestID("cellxgene-logo"));
  });

  it("renders the the Privacy Page", async () => {
    await goToPage(`${TEST_URL}${ROUTES.PRIVACY}`);
    page.screenshot({ path: `./tmp-screenshots/privacy-${Date.now()}.png` });

    await expect(page).toHaveSelector(getText("Privacy Policy"));
    await expect(page).toHaveSelector(getTestID("cellxgene-logo"));
  });
});
