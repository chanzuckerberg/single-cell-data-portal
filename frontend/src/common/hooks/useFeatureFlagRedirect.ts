// Core dependencies
import { useRouter } from "next/router";
import { useEffect, useState } from "react";
// App dependencies
import { ROUTES } from "src/common/constants/routes";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";

/**
 * Redirect to given route if feature is not available to user.
 * @param featureFlag - Name of feature to check if available to user.
 * @param redirectRoute - Route to redirect user to if feature is not available to them.
 * @returns True if feature is available to user.
 */
export function useFeatureFlagRedirect(
  featureFlag: FEATURES,
  redirectRoute: ROUTES
): boolean {
  const router = useRouter();
  const { isReady } = router;

  /* Flag indicating if feature is available to user. */
  const [isFilterEnabled, setIsFilterEnabled] = useState<boolean>(false);

  /* Update state of enabled flag and redirect user if feature is not available to them. */
  useEffect(() => {
    const enabled = get(FEATURES.FILTER) === BOOLEAN.TRUE;
    setIsFilterEnabled(enabled);
    if (!enabled) {
      router.push(ROUTES.HOMEPAGE);
    }
  }, [featureFlag, isFilterEnabled, isReady, redirectRoute, router]);

  return isFilterEnabled;
}
