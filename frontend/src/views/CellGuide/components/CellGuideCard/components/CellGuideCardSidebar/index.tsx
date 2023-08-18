import { SearchBarWrapper } from "./style";
import CellGuideCardSearchBar from "../../../CellGuideCardSearchBar";
import {
  CellGuideSidebarWrapper,
  StickyWrapper,
  StickySidebarStyle,
} from "./style";
import { NavigationJumpTo } from "@czi-sds/components";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import { Global } from "@emotion/react";
import { useEffect, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";

export const CELL_GUIDE_CARD_NAVIGATION_SIDEBAR =
  "cell-guide-card-navigation-sidebar";

export default function CellGuideCardSidebar({
  items,
}: {
  items: { elementRef: React.MutableRefObject<null>; title: string }[];
}): JSX.Element {
  /**
   * (thuang): SDS `NavigationJumpTo` currently assumes `window` is always available,
   * which is not the case for SSR. We need to use `isClient` to prevent the component
   * from rendering on the server.
   *
   * Remove this after SDS fixes the bug
   * PR: https://github.com/chanzuckerberg/sci-components/pull/514
   */
  const [isClient, setIsClient] = useState(false);

  useEffect(() => {
    setIsClient(true);
  }, []);

  return (
    <CellGuideSidebarWrapper data-testid={CELL_GUIDE_CARD_NAVIGATION_SIDEBAR}>
      <Global styles={StickySidebarStyle} />
      <StickyWrapper>
        <SearchBarWrapper>
          <CellGuideCardSearchBar />
        </SearchBarWrapper>
        {isClient && (
          <NavigationJumpTo
            items={items}
            offsetTop={HEADER_HEIGHT_PX} // Used to account for header height
            onClick={(event) => {
              track(EVENTS.CG_VIEW_SECTION, {
                section: (event.target as HTMLElement).innerText,
              });
            }}
          />
        )}
      </StickyWrapper>
    </CellGuideSidebarWrapper>
  );
}
