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

export default function CellCardSidebar({
  items,
}: {
  items: { elementRef: React.MutableRefObject<null>; title: string }[];
}): JSX.Element {
  return (
    <CellCardsSidebarWrapper>
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
