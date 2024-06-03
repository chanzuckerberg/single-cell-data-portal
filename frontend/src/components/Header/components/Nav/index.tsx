import React from "react";
import Link from "next/link";
import { ROUTES } from "src/common/constants/routes";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import NavDivider from "src/components/Header/components/Nav/components/NavDivider";
import { useConnect } from "src/components/Header/connect";
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
  const { isRouteActive } = useConnect();
  return (
    <>
      <NavWrapper className={className}>
        <NavSection>
          <NavSectionTitle>Application</NavSectionTitle>
          <NavItemContainer>
            <LinkWrapper
              isActive={isRouteActive(pathname, ROUTES.COLLECTIONS)}
              data-testid="collections-link"
              onClick={() => {
                track(EVENTS.COLLECTIONS_CLICK_NAV);
              }}
            >
              <Link href={ROUTES.COLLECTIONS}>Collections</Link>
            </LinkWrapper>
            <LinkWrapper
              isActive={isRouteActive(pathname, ROUTES.DATASETS)}
              onClick={() => {
                track(EVENTS.DATASETS_CLICK_NAV);
              }}
            >
              <Link href={ROUTES.DATASETS}>Datasets</Link>
            </LinkWrapper>
            <LinkWrapper
              isActive={isRouteActive(pathname, ROUTES.WHERE_IS_MY_GENE)}
            >
              <Link
                href={ROUTES.WHERE_IS_MY_GENE}
                onClick={() => {
                  track(EVENTS.WMG_CLICK_NAV);
                }}
              >
                Gene Expression
              </Link>
            </LinkWrapper>
            <LinkWrapper
              isActive={isRouteActive(pathname, ROUTES.CELL_GUIDE)}
              onClick={() => {
                track(EVENTS.CELL_GUIDE_CLICK_NAV);
              }}
            >
              <Link href={ROUTES.CELL_GUIDE}>Cell Guide</Link>
              <BetaChip label="Beta" size="small" />
            </LinkWrapper>
            {/* Uncomment this when DE is ready for release */}
            {/* <LinkWrapper
              isActive={isRouteActive(pathname, ROUTES.DE)}
              onClick={() => {
                track(EVENTS.DE_CLICK_NAV);
              }}
            >
              <Link href={ROUTES.DE}>Differential Expression</Link>
              <BetaChip label="Beta" size="small" />
            </LinkWrapper> */}
          </NavItemContainer>
        </NavSection>
        <NavDivider />
        <NavSection>
          <NavSectionTitle>Census</NavSectionTitle>
          <NavItemContainer>
            <LinkWrapper isActive={false}>
              <a
                onClick={() => track(EVENTS.CENSUS_DOCUMENTATION_CLICK_NAV)}
                href={CENSUS_LINK}
                rel="noopener noreferrer"
              >
                API
              </a>
            </LinkWrapper>
            <LinkWrapper
              isActive={isRouteActive(pathname, ROUTES.CENSUS_DIRECTORY)}
              onClick={() => {
                track(EVENTS.CENSUS_DIRECTORY_CLICK_NAV);
              }}
            >
              <Link href={ROUTES.CENSUS_DIRECTORY}>Models</Link>
            </LinkWrapper>
          </NavItemContainer>
        </NavSection>
      </NavWrapper>
    </>
  );
}
