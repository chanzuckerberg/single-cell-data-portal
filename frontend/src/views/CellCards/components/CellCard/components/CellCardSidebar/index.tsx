import { SearchBarWrapper } from "./style";
import CellCardSearchBar from "../../../CellCardSearchBar";
import { CellCardsSidebarWrapper, StickyWrapper } from "./style";
import { NavigationJumpTo } from "@czi-sds/components";

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
          <NavigationJumpTo
            items={items}
            content={undefined}
            rel={undefined}
            rev={undefined}
            slotProps={undefined}
            slots={undefined}
          />
        </StickyWrapper>
      </CellCardsSidebarWrapper>
    </>
  );
}
