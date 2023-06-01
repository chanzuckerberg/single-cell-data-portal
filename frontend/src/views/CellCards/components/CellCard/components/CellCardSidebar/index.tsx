import { SearchBarWrapper } from "./style";
import CellCardSearchBar from "../../../CellCardSearchBar";
import {
  CellCardsSidebarWrapper,
  StickyWrapper,
  StickySidebarStyle,
} from "./style";
import { NavigationJumpTo } from "@czi-sds/components";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import { Global } from "@emotion/react";

export const CELL_CARD_NAVIGATION_SIDEBAR = "cell-card-navigation-sidebar";

export default function CellCardSidebar({
  items,
}: {
  items: { elementRef: React.MutableRefObject<null>; title: string }[];
}): JSX.Element {
  return (
    <CellCardsSidebarWrapper data-testid={CELL_CARD_NAVIGATION_SIDEBAR}>
      <Global styles={StickySidebarStyle} />
      <StickyWrapper>
        <SearchBarWrapper>
          <CellCardSearchBar />
        </SearchBarWrapper>
        <NavigationJumpTo
          items={items}
          offsetTop={HEADER_HEIGHT_PX} // Used to account for header height
        />
      </StickyWrapper>
    </CellCardsSidebarWrapper>
  );
}
