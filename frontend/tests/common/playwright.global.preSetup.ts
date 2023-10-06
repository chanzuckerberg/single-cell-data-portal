import { test } from "@playwright/test";
import { getAccessToken } from "tests/utils/helpers";

const { describe } = test;

/**
 * (thuang): This file sets up access token needed to get through gated rdev stacks
 * as needed
 */
describe("global preSetup", () => {
  test("Get access token", async ({ request }) => {
    await getAccessToken(request);
  });
});
