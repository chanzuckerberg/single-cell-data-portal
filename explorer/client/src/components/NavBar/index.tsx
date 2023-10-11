import React from "react";
import { AnchorButton, Menu, MenuItem, Position } from "@blueprintjs/core";
import { Popover2 } from "@blueprintjs/popover2";
import { IconNames as CXGIconNames } from "../icon";
import { track } from "../../analytics";
import { EVENTS } from "../../analytics/events";
import { ROUTES } from "./routes";
import Icon from "../icon/icon";
import {
  BetaChip,
  Left,
  LinkWrapper,
  MainWrapper,
  Nav,
  Right,
  Wrapper,
} from "./style";
import NavDivider from "./components/NavDivider";

function handleMenuClick() {
  track(EVENTS.EXPLORER_MENU_BUTTON_CLICKED);
}

const CENSUS_LINK = "https://cellxgene-census.readthedocs.io/en/latest";

interface HeaderProps {
  tosURL?: string;
  privacyURL?: string;
}
const Header = (props: HeaderProps) => {
  const { tosURL, privacyURL } = props;
  return (
    <Wrapper data-testid="header">
      <MainWrapper>
        <Left>
          <a href={ROUTES.HOMEPAGE}>
            <svg
              width="23"
              height="23"
              viewBox="0 0 23 23"
              fill="none"
              xmlns="http://www.w3.org/2000/svg"
            >
              <path
                d="M3.10758 5.97949H0.521711C0.233578 5.97949 0 6.21307 0 6.5012V22.1616C0 22.4497 0.233578 22.6833 0.521711 22.6833H3.10758C3.39571 22.6833 3.62929 22.4497 3.62929 22.1616V6.5012C3.62929 6.21307 3.39571 5.97949 3.10758 5.97949Z"
                fill="white"
              />
              <path
                d="M22.6741 22.1613C22.6741 22.4471 22.4382 22.683 22.1523 22.683H6.50102C6.21521 22.683 5.97931 22.4471 5.97931 22.1613V6.50999C5.97931 6.22419 6.21521 5.98828 6.50102 5.98828H22.1523C22.4382 5.98828 22.6741 6.22419 22.6741 6.50999V22.1613ZM9.6086 18.532C9.6086 18.8178 9.84451 19.0537 10.1303 19.0537H18.5231C18.8089 19.0537 19.0448 18.8178 19.0448 18.532V10.1393C19.0448 9.85348 18.8089 9.61757 18.5231 9.61757H10.1303C9.84451 9.61757 9.6086 9.85348 9.6086 10.1393V18.532Z"
                fill="white"
              />
              <path
                d="M22.1569 0H6.4965C6.20837 0 5.97479 0.233578 5.97479 0.521711V3.10758C5.97479 3.39571 6.20837 3.62929 6.4965 3.62929H22.1569C22.445 3.62929 22.6786 3.39571 22.6786 3.10758V0.521711C22.6786 0.233578 22.445 0 22.1569 0Z"
                fill="#8282FF"
              />
            </svg>
          </a>
          <Nav>
            <LinkWrapper>
              <AnchorButton
                active={false}
                href={ROUTES.COLLECTIONS}
                minimal
                text="Collections"
                onClick={handleCollectionsClick}
              />
            </LinkWrapper>
            <LinkWrapper>
              <AnchorButton
                active={false}
                href={ROUTES.DATASETS}
                minimal
                text="Datasets"
                onClick={handleDatasetsClick}
              />
            </LinkWrapper>
            <LinkWrapper>
              <AnchorButton
                active={false}
                href={ROUTES.WHERE_IS_MY_GENE}
                minimal
                text="Gene Expression"
                onClick={handleWMGClick}
              />
            </LinkWrapper>
            <LinkWrapper>
              <AnchorButton
                active={false}
                href={ROUTES.CELL_GUIDE}
                minimal
                text="Cell Guide"
                onClick={handleCellGuideClick}
              />
              <BetaChip label="Beta" size="small" />
            </LinkWrapper>
            <NavDivider />
            <LinkWrapper>
              <AnchorButton
                href={CENSUS_LINK}
                minimal
                onClick={handleCensusClick}
                rel="noopener"
                target="_self"
                text="Census"
              />
            </LinkWrapper>
          </Nav>
        </Left>
        <Right>
          <LinkWrapper>
            <Popover2
              content={
                <Menu>
                  <MenuItem
                    href={ROUTES.DOCS}
                    target="_blank"
                    text="Documentation"
                    rel="noopener"
                    onClick={handleDocumentationClick}
                  />
                  <MenuItem
                    href="https://join-cellxgene-users.herokuapp.com/"
                    icon={<Icon icon={CXGIconNames.SLACK} />}
                    target="_blank"
                    text="Chat"
                    rel="noopener"
                  />
                  {(tosURL || privacyURL) && (
                    <MenuItem
                      icon={<Icon icon={CXGIconNames.ABOUT} />}
                      popoverProps={{ openOnTargetFocus: false }}
                      text="About cellxgene"
                    >
                      {tosURL && (
                        <MenuItem
                          href={tosURL}
                          rel="noopener"
                          target="_blank"
                          text="Terms of Service"
                        />
                      )}
                      {privacyURL && (
                        <MenuItem
                          href={privacyURL}
                          rel="noopener"
                          target="_blank"
                          text="Privacy Policy"
                        />
                      )}
                    </MenuItem>
                  )}
                </Menu>
              }
              position={Position.BOTTOM_LEFT}
              modifiers={{
                preventOverflow: { enabled: false },
                hide: { enabled: false },
              }}
            >
              <AnchorButton
                active={false}
                data-testid="menu"
                minimal
                text="Help & Documentation"
                rightIcon="chevron-down"
                onClick={handleMenuClick}
              />
            </Popover2>
          </LinkWrapper>
        </Right>
      </MainWrapper>
    </Wrapper>
  );

  function handleWMGClick(): void {
    track(EVENTS.WMG_CLICK_NAV);
  }
  function handleCellGuideClick(): void {
    track(EVENTS.CELL_GUIDE_CLICK_NAV);
  }
  function handleDatasetsClick(): void {
    track(EVENTS.DATASETS_CLICK_NAV);
  }
  function handleCollectionsClick(): void {
    track(EVENTS.COLLECTIONS_CLICK_NAV);
  }
  function handleCensusClick(): void {
    track(EVENTS.CENSUS_CLICK_NAV);
  }
  function handleDocumentationClick(): void {
    track(EVENTS.DOCUMENTATION_CLICK_NAV);
  }
};
Header.defaultProps = {
  tosURL: undefined,
  privacyURL: undefined,
};
export default Header;
