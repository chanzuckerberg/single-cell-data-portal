// Core dependencies
import { useRouter } from "next/router";
import { useEffect, useState } from "react";
// App dependencies
import { ROUTES } from "src/common/constants/routes";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";

/**
 * Determine if feature is available to user. Redirect to route, if specified.
 * @param featureFlag - Name of feature to check if available to user.
 * @param redirectRoute - Route to redirect user to - if any - if feature is not available to them.
 * @returns True if feature is available to user.
 */
export function useFeatureFlag(
  featureFlag: FEATURES,
  redirectRoute?: ROUTES
): boolean {
  const router = useRouter();

  /* Flag indicating if feature is available to user. */
  const [isEnabled, setIsEnabled] = useState<boolean>(false);

  /* Update state of enabled flag and redirect user if feature is not available to them. */
  useEffect(() => {
    // Force enable of filter - start. Remove with #1718.
    if (featureFlag === FEATURES.FILTER) {
      setIsEnabled(true);
      return;
    }
    // Force enable of filter - end.
    const enabled = get(featureFlag) === BOOLEAN.TRUE;
    setIsEnabled(enabled);
    if (!enabled && redirectRoute) {
      router.push(redirectRoute);
    }
  }, [featureFlag, redirectRoute, router]);

  return isEnabled;
}
