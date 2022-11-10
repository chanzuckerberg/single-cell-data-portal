import { TEST_URL } from "./constants";

const featureFlags = {
  origins: [
    {
      localStorage: [
        { name: "cxg-ff-gs", value: "yes" },
        { name: "cxg-ff-curator", value: "yes" },
        { name: "cxg-ff-rc", value: "yes" },
      ],
      origin: TEST_URL,
    },
  ],
};

export default featureFlags;
