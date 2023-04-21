import { AnchorButton } from "@blueprintjs/core";
import Link from "next/link";
import { useRouter } from "next/router";
import { FC, useRef, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ROUTES } from "src/common/constants/routes";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useUserInfo } from "src/common/queries/auth";
import AuthButtons from "src/components/Header/components/AuthButtons";
import { HomepageLink } from "../common/HomepageLink";
import {
  DesktopHomeLink,
  HiringLink,
  Left,
  LinkWrapper,
  MainWrapper,
  MobileHomeLink,
  MobileMenuButton,
  MobileMenuButtonBar,
  MobileNavTray,
  MobileNavWrapper,
  Nav,
  Right,
  Wrapper,
} from "./style";

const LandingHeader: FC = () => {
  const isCurator = get(FEATURES.CURATOR) === BOOLEAN.TRUE;
  const { data: userInfo } = useUserInfo(isCurator);
  const { pathname } = useRouter();
  const isMyCollectionsShown = userInfo?.name && isCurator;

  const mobileNavTray = useRef<HTMLDivElement>(null);
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false);

  function mobileNavHandler(mobileMenuOpen: boolean) {
    if (!mobileMenuOpen) {
      setMobileMenuOpen(true);
      document.documentElement.style.overflowY = "hidden";
    } else {
      document.documentElement.style.overflowY = "visible";
      setMobileMenuOpen(false);
    }
  }

  return (
    <MobileNavWrapper>
      <MobileHomeLink>
        <HomepageLink />
      </MobileHomeLink>
      <MobileMenuButton onClick={() => mobileNavHandler(mobileMenuOpen)}>
        <MobileMenuButtonBar className={mobileMenuOpen ? "open" : ""} />
        <MobileMenuButtonBar className={mobileMenuOpen ? "open" : ""} />
        <MobileMenuButtonBar className={mobileMenuOpen ? "open" : ""} />
      </MobileMenuButton>
      <MobileNavTray
        className={`${mobileMenuOpen ? "active" : ""}`}
        ref={mobileNavTray}
      >
        <Wrapper>
          <MainWrapper>
            <Left>
              <DesktopHomeLink>
                <HomepageLink />
              </DesktopHomeLink>
              <Nav>
                <LinkWrapper>
                  <Link href={ROUTES.COLLECTIONS} passHref>
                    <AnchorButton
                      data-testid="collection-link"
                      onClick={() => {
                        track(EVENTS.COLLECTIONS_CLICK_NAV);
                      }}
                      active={isRouteActive(pathname, ROUTES.COLLECTIONS)}
                      href="passHref"
                      minimal
                      text="Collections"
                    />
                  </Link>
                </LinkWrapper>
                <LinkWrapper>
                  <Link href={ROUTES.DATASETS} passHref>
                    <AnchorButton
                      onClick={() => {
                        track(EVENTS.DATASETS_CLICK_NAV);
                      }}
                      active={isRouteActive(pathname, ROUTES.DATASETS)}
                      href="passHref"
                      minimal
                      text="Datasets"
                    />
                  </Link>
                </LinkWrapper>
                <LinkWrapper>
                  <Link href={ROUTES.WHERE_IS_MY_GENE} passHref>
                    <AnchorButton
                      active={isRouteActive(pathname, ROUTES.WHERE_IS_MY_GENE)}
                      href="passHref"
                      minimal
                      text="Gene Expression"
                      onClick={handleWMGClick}
                    />
                  </Link>
                </LinkWrapper>
                <LinkWrapper>
                  <Link href={ROUTES.CELL_CARDS} passHref>
                    <AnchorButton
                      active={isRouteActive(pathname, ROUTES.CELL_CARDS)}
                      href="passHref"
                      minimal
                      text="Cell Cards"
                      onClick={handleCellCardsClick}
                    />
                  </Link>
                </LinkWrapper>
              </Nav>
            </Left>
            <Right>
              {pathname === ROUTES.HOMEPAGE && (
                <HiringLink
                  onClick={() => {
                    track(EVENTS.BROWSE_CAREERS_CLICKED, {
                      button: "we're hiring",
                    });
                  }}
                  href="https://chanzuckerberg.com/careers/career-opportunities/?team=data,design,engineering,product,technical-program-management&initiative=science&gh_src=20d9f28d1us"
                  target="_blank"
                  rel="noopener"
                >
                  We&apos;re Hiring!
                </HiringLink>
              )}
              {isMyCollectionsShown && (
                <LinkWrapper>
                  <Link href={ROUTES.MY_COLLECTIONS} passHref>
                    <AnchorButton
                      active={isRouteActive(pathname, ROUTES.MY_COLLECTIONS)}
                      href="passHref"
                      minimal
                      text="My Collections"
                    />
                  </Link>
                </LinkWrapper>
              )}
              <LinkWrapper>
                <AnchorButton
                  onClick={() => {
                    track(EVENTS.DOCUMENTATION_CLICK_NAV);
                  }}
                  active={isRouteActive(pathname, ROUTES.DOCS)}
                  href={ROUTES.DOCS}
                  rel="noopener"
                  target="_blank"
                  minimal
                  text="Help & Documentation"
                />
              </LinkWrapper>
              <AuthButtons />
            </Right>
          </MainWrapper>
        </Wrapper>
      </MobileNavTray>
    </MobileNavWrapper>
  );

  function handleWMGClick() {
    track(EVENTS.WMG_CLICK_NAV);
  }

  function handleCellCardsClick() {
    track(EVENTS.CELL_CARDS_CLICK_NAV);
  }
};

/**
 * Returns true if current path is equal to the route.
 * @param path
 * @param route
 * @returns true if the current path is the route
 */
function isRouteActive(path: string, route: ROUTES): boolean {
  return path === route;
}

export default LandingHeader;
