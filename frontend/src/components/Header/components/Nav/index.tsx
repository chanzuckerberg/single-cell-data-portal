import React from "react";
import Link from "next/link";
import { ROUTES } from "src/common/constants/routes";
import { AnchorButton } from "@blueprintjs/core";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import NavDivider from "src/components/Header/components/Nav/components/NavDivider";
import { isRouteActive } from "src/components/Header";
import {
  LinkWrapper,
  NavItemContainer,
  NavSection,
  NavSectionTitle,
  Nav as NavWrapper,
} from "./style";
import { BetaChip } from "../../style";
import { CENSUS_LINK } from "./constants";
import { Props } from "./types";

export default function Nav({ className, pathname }: Props): JSX.Element {
  return (
    <>
      <NavWrapper className={className}>
        <NavSection>
          <NavSectionTitle>Application</NavSectionTitle>
          <NavItemContainer>
            <LinkWrapper>
              <Link href={ROUTES.COLLECTIONS} passHref>
                <AnchorButton
                  active={isRouteActive(pathname, ROUTES.COLLECTIONS)}
                  data-testid="collections-link"
                  href="passHref"
                  minimal
                  onClick={() => {
                    track(EVENTS.COLLECTIONS_CLICK_NAV);
                  }}
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
                  onClick={() => {
                    track(EVENTS.DATASETS_CLICK_NAV);
                  }}
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
                  onClick={() => {
                    track(EVENTS.WMG_CLICK_NAV);
                  }}
                  text="Gene Expression"
                />
              </Link>
            </LinkWrapper>
            <LinkWrapper>
              <Link href={ROUTES.CELL_GUIDE} passHref>
                <AnchorButton
                  active={isRouteActive(pathname, ROUTES.CELL_GUIDE)}
                  href="passHref"
                  minimal
                  onClick={() => {
                    track(EVENTS.CELL_GUIDE_CLICK_NAV);
                  }}
                  text="Cell Guide"
                />
              </Link>
              <BetaChip label="Beta" size="small" />
            </LinkWrapper>
          </NavItemContainer>
        </NavSection>
        <NavDivider />
        <NavSection>
          <NavSectionTitle>Census</NavSectionTitle>
          <NavItemContainer>
            <LinkWrapper>
              <AnchorButton
                href={CENSUS_LINK}
                minimal
                onClick={() => {
                  track(EVENTS.CENSUS_DOCUMENTATION_CLICK_NAV);
                }}
                rel="noopener"
                target="_self"
                text="API"
              />
            </LinkWrapper>
            <LinkWrapper>
              <Link href={ROUTES.CENSUS_DIRECTORY} passHref>
                <AnchorButton
                  active={isRouteActive(pathname, ROUTES.CENSUS_DIRECTORY)}
                  href="passHref"
                  minimal
                  onClick={() => {
                    track(EVENTS.CENSUS_DIRECTORY_CLICK_NAV);
                  }}
                  text="Models"
                />
              </Link>
            </LinkWrapper>
          </NavItemContainer>
        </NavSection>
      </NavWrapper>
    </>
  );
}
