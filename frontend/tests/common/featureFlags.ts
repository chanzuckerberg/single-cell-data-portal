import { TEST_URL } from "./constants";

/**
 * (thuang): Set curator flag to true to enable login.
 */
export function getFeatureFlags({
  geneSet = true,
  curator = false,
}: { geneSet?: boolean; curator?: boolean } = {}): {
  cookies: Array<{
    name: string;

    value: string;

    /**
     * domain and path are required
     */
    domain: string;

    /**
     * domain and path are required
     */
    path: string;

    /**
     * Unix time in seconds.
     */
    expires: number;

    httpOnly: boolean;

    secure: boolean;

    /**
     * sameSite flag
     */
    sameSite: "Strict" | "Lax" | "None";
  }>;
  origins: Array<{
    origin: string;

    localStorage: Array<{
      name: string;

      value: string;
    }>;
  }>;
} {
  return {
    cookies: [],
    origins: [
      {
        localStorage: [
          { name: "cxg-ff-gs", value: geneSet ? "yes" : "no" },
          { name: "cxg-ff-curator", value: curator ? "yes" : "no" },
        ],
        origin: TEST_URL,
      },
    ],
  };
}
