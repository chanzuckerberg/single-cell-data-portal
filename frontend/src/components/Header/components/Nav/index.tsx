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

export const CENSUS_LINK = "https://chanzuckerberg.github.io/cellxgene-census/";

interface Props {
  className?: string;
  pathname: string;
}

export default function Nav({ className, pathname }: Props): JSX.Element {
  const isCensusDirectory = isRouteActive(pathname, ROUTES.CENSUS_DIRECTORY);

  return !isCensusDirectory ? (
    <NavWrapper className={className}>
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
      <NavDivider />
      <LinkWrapper>
        <AnchorButton
          href={CENSUS_LINK}
          minimal
          onClick={() => {
            track(EVENTS.CENSUS_CLICK_NAV);
          }}
          rel="noopener"
          target="_self"
          text="Census"
        />
      </LinkWrapper>
    </NavWrapper>
  ) : (
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
                  track(EVENTS.CENSUS_CLICK_NAV);
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
                    // TODO ANALYTICS
                  }}
                  text="Spotlight"
                />
              </Link>
            </LinkWrapper>
          </NavItemContainer>
        </NavSection>
      </NavWrapper>
    </>
  );
}
