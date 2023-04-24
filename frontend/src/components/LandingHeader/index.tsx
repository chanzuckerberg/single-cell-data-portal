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
  Left,
  MainWrapper,
  MobileHomeLink,
  MobileMenuButton,
  MobileMenuButtonBar,
  MobileNavTray,
  MobileNavWrapper,
  Navigation as Nav,
  Right,
  Wrapper,
} from "./style";
import { LinkWrapper } from "src/components/Header/components/Nav/style";

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
              <Nav pathname={pathname} />
            </Left>
            <Right>
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
