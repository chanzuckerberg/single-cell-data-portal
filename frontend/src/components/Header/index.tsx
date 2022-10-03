import { AnchorButton } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { Popper } from "@material-ui/core";
import { DefaultMenuSelectOption, InputDropdown, MenuSelect } from "czifui";
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
  `;

  return (
    <Wrapper>
      <MainWrapper>
        <Left>
          <HomepageLink />
          <Nav>
            <LinkWrapper>
              <Link href={ROUTES.COLLECTIONS} passHref>
                <AnchorButton
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

          <InputDropdown
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
    let link = newValue!.value;

    if (newValue!.id === 1) {
      track(EVENTS.DOCUMENTATION_CLICK_NAV);

      link = isRouteActive(pathname, ROUTES.WHERE_IS_MY_GENE)
        ? ROUTES.WMG_DOCS
        : ROUTES.PUBLISHED_DATA_DOCS;
    }

    window.open(link, "_blank", "noopener,noreferrer");

    setDropdownOpen(false);
  }

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
