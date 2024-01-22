import { DefaultMenuSelectOption } from "@czi-sds/components";
import { useRouter } from "next/router";
import { useRef, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ROUTES } from "src/common/constants/routes";

export function useConnect() {
  /**
   * Retrieves the current pathname from the router.
   * @returns The current pathname.
   */
  const { pathname } = useRouter();

  /**
   * A reference to the dropdown element.
   */
  const dropdownRef = useRef(null);

  /**
   * The anchor element for the dropdown menu.
   */
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(
    dropdownRef.current
  );

  /**
   * Represents the state of the dropdown menu.
   */
  const [dropdownOpen, setDropdownOpen] = useState(false);

  /**
   * Handles the opening of the help dropdown menu.
   *
   * @param event - The mouse event that triggered the function.
   */
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

  /**
   * Handles the click event for the help dropdown menu.
   * @param {React.ChangeEvent<unknown>} _ - The event object.
   * @param {DefaultMenuSelectOption & { id: number; value: string } | null} newValue - The selected option from the dropdown menu.
   */
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

  /**
   * Closes the help dropdown.
   */
  function handleHelpClose() {
    setDropdownOpen(false);
  }

  /**
   * Returns true if current path is equal to the route.
   * @param path
   * @param route
   * @returns true if the current path is the route
   */
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
