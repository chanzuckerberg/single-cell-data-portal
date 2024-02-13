import { DefaultMenuSelectOption } from "@czi-sds/components";
import { useRouter } from "next/router";
import { useRef, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ROUTES } from "src/common/constants/routes";

export function useConnect() {
  const { pathname } = useRouter();

  const dropdownRef = useRef(null);

  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(
    dropdownRef.current
  );

  const [dropdownOpen, setDropdownOpen] = useState(false);

  function handleHelpOpen(event: React.MouseEvent<HTMLElement>) {
    if (!anchorEl) {
      setAnchorEl(event.currentTarget);
    }

    if (dropdownOpen) {
      setDropdownOpen(false);
    } else {
      setDropdownOpen(true);
    }
  }

  function handleHelpClick(
    _: React.ChangeEvent<unknown>,
    newValue: (DefaultMenuSelectOption & { id: number; value: string }) | null
  ) {
    if (!newValue) return;
    let link = newValue.value;

    if (newValue.id === 1) {
      track(EVENTS.DOCUMENTATION_CLICK_NAV);

      link = isRouteActive(pathname, ROUTES.WHERE_IS_MY_GENE)
        ? ROUTES.WMG_DOCS
        : ROUTES.PUBLISHED_DATA_DOCS;
    }

    window.open(link, "_blank", "noopener,noreferrer");
    setDropdownOpen(false);
  }

  function handleHelpClose() {
    setDropdownOpen(false);
  }

  function isRouteActive(path: string, route: ROUTES): boolean {
    return path === route;
  }

  return {
    pathname,
    anchorEl,
    dropdownOpen,
    dropdownRef,
    handleHelpOpen,
    handleHelpClick,
    handleHelpClose,
    isRouteActive,
  };
}
