import { AnchorButton } from "@blueprintjs/core";
import Link from "next/link";
import { useRouter } from "next/router";
import { FC } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ROUTES } from "src/common/constants/routes";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useUserInfo } from "src/common/queries/auth";
import { HomepageLink } from "../common/HomepageLink";
import AuthButtons from "./components/AuthButtons";
import {
  BetaChip,
  Left,
  LinkWrapper,
  MainWrapper,
  Nav,
  Right,
  Wrapper,
} from "./style";

const Header: FC = () => {
  const isCurator = get(FEATURES.CURATOR) === BOOLEAN.TRUE;
  const { data: userInfo } = useUserInfo(isCurator);
  const { pathname } = useRouter();
  const isMyCollectionsShown = userInfo?.name && isCurator;

  return (
    <Wrapper>
      <MainWrapper>
        <Left>
          <HomepageLink />
          <Nav>
            <LinkWrapper>
              <Link href={ROUTES.COLLECTIONS} passHref>
                <AnchorButton
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
                  text="scExpression"
                  onClick={handleWMGClick}
                />
              </Link>
              <BetaChip label="Beta" size="small" />
            </LinkWrapper>
          </Nav>
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
  );

  function handleWMGClick() {
    track(EVENTS.WMG_CLICK_NAV);
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

export default Header;
