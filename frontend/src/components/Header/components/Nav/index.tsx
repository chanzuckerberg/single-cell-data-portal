import React from "react";
import Link from "next/link";
import { ROUTES } from "src/common/constants/routes";
import { AnchorButton } from "@blueprintjs/core";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import NavDivider from "src/components/Header/components/Nav/components/NavDivider";
import { isRouteActive } from "src/components/Header";
import { LinkWrapper, Nav as NavWrapper } from "./style";

const CENSUS_LINK = "https://cellxgene-census.readthedocs.io/en/latest";

interface Props {
  className?: string;
  pathname: string;
}

export default function Nav({ className, pathname }: Props): JSX.Element {
  return (
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
      {/* <LinkWrapper>
        <Link href={ROUTES.CELL_CARDS} passHref>
          <AnchorButton
            active={isRouteActive(pathname, ROUTES.CELL_CARDS)}
            href="passHref"
            minimal
            onClick={() => {
              track(EVENTS.CELL_CARDS_CLICK_NAV);
            }}
            text="Cell Cards"
          />
        </Link>
      </LinkWrapper> */}
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
  );
}
