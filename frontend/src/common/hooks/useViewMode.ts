import { useMemo } from "react";
import { get } from "src/common/featureFlags";
import { BOOLEAN } from "src/common/localStorage/set";
import { FEATURES } from "src/common/featureFlags/features";
import { useUserInfo } from "src/common/queries/auth";
import { QueryStatus } from "react-query";

/**
 * View mode (view is either default or authenticated / curator mode).
 */
export enum VIEW_MODE {
  CURATOR = "CURATOR",
  DEFAULT = "DEFAULT",
}

interface UseViewMode {
  mode: VIEW_MODE;
  status: QueryStatus;
}

/**
 * Determines the view mode and user info query status for the index view of collections and datasets.
 * The view mode can be either "default" or "curator" mode and is used to enable or disable features in the
 * collections view, such as the collections list and filtering features, as well as toggle the collections and
 * datasets query function.
 * The user info query status is responsible for enabling the collections and datasets query,
 * both of which have an infinite stale time. It is important to fetch these queries only when the user info query is
 * complete and the mode is determined.
 * @returns view mode and status of user info query.
 */
export function useViewMode(): UseViewMode {
  const isCuratorEnabled = get(FEATURES.CURATOR) === BOOLEAN.TRUE;
  const { data: userInfo, status } = useUserInfo(isCuratorEnabled);
  const isAuth = !!userInfo?.name && status === "success";

  // Determine view mode.
  const mode = useMemo((): VIEW_MODE => {
    return isAuth ? VIEW_MODE.CURATOR : VIEW_MODE.DEFAULT;
  }, [isAuth]);

  return { mode, status };
}
