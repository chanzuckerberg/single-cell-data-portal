import { SearchBarWrapper } from "./style";
import CellCardSearchBar from "../../../CellCardSearchBar";
import { CellCardsSidebarWrapper, StickyWrapper } from "./style";
import { NavigationJumpTo } from "@czi-sds/components";
import { isSSR } from "src/common/utils/isSSR";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

export default function CellCardSidebar({
  items,
}: {
  items: { elementRef: React.MutableRefObject<null>; title: string }[];
}): JSX.Element {
  return (
    <>
      <style>
        {`
          /* Hack because main has a global overflow CSS prop which interferes with sticky sidebar */
          main {
            overflow: unset !important;
          }
        `}
      </style>
      <CellCardsSidebarWrapper>
        <StickyWrapper>
          <SearchBarWrapper>
            <CellCardSearchBar />
          </SearchBarWrapper>

          {/* Unsure why some of these props are required as they aren't needed, SDS might fix in future a future version */}
          {/* Also there is a SDS bug where you get "window is not defined" error without checking isSSR */}
          {!isSSR() && (
            <NavigationJumpTo
              items={items}
              offsetTop={HEADER_HEIGHT_PX} // Used to account for header height
              content={undefined}
              rel={undefined}
              rev={undefined}
              slotProps={undefined}
              slots={undefined}
            />
          )}
        </StickyWrapper>
      </CellCardsSidebarWrapper>
    </>
  );
}
