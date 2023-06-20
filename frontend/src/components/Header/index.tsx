import styled from "@emotion/styled";
import { Popper } from "@mui/material";
import { DefaultMenuSelectOption, MenuSelect } from "@czi-sds/components";
import { useRouter } from "next/router";
import { FC, useRef, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ROUTES } from "src/common/constants/routes";
// import { noop } from "src/common/constants/utils";
import { HomepageLink } from "../common/HomepageLink";
import AuthButtons from "./components/AuthButtons";
import {
  Left,
  MainWrapper,
  Right,
  StyledInputDropdown,
  Wrapper,
} from "./style";
import Nav from "src/components/Header/components/Nav";

const Header: FC = () => {
  const { pathname } = useRouter();

  const dropdownRef = useRef(null);

  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(
    dropdownRef.current
  );
  const [dropdownOpen, setDropdownOpen] = useState(false);

  const FOOTER_OPTIONS = [
    { id: 1, name: "Documentation", value: ROUTES.WMG_DOCS },
    { id: 2, name: "Contact", value: "mailto:cellxgene@chanzuckerberg.com" },
    { id: 3, name: "Terms of Service", value: ROUTES.TOS },
    { id: 4, name: "Privacy Policy", value: ROUTES.PRIVACY },
  ];

  const StyledPopper = styled(Popper)`
    height: 300px;
    z-index: 99;
    margin-top: -26px !important;
  `;

  return (
    <Wrapper>
      <MainWrapper>
        <Left>
          <HomepageLink />
          <Nav pathname={pathname} />
        </Left>
        <Right>
          <StyledInputDropdown
            disabled={false}
            label="Help & Documentation"
            sdsStage="default"
            onClick={handleHelpOpen}
            sdsStyle="minimal"
            sdsType="singleSelect"
            data-testid="InputDropdown"
          />

          <StyledPopper
            ref={dropdownRef}
            open={dropdownOpen}
            anchorEl={anchorEl}
          >
            <MenuSelect
              search={false}
              options={FOOTER_OPTIONS}
              onChange={handleHelpClick}
              onClose={handleHelpClose}
            />
          </StyledPopper>

          <AuthButtons />
        </Right>
      </MainWrapper>
    </Wrapper>
  );

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
};

/**
 * Returns true if current path is equal to the route.
 * @param path
 * @param route
 * @returns true if the current path is the route
 */
export function isRouteActive(path: string, route: ROUTES): boolean {
  return path === route;
}

export default Header;
