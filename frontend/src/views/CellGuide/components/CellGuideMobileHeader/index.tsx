import { useContext, useState } from "react";
import CellGuideCardSearchBar from "../CellGuideCardSearchBar";
import { ButtonIcon } from "@czi-sds/components";
import {
  MobileHeader,
  MobileHeaderWrapper,
  MobilePageNavWrapper,
  MobileSearchBarWrapper,
  StyledTitle,
} from "./style";
import { DispatchContext, StateContext } from "../../common/store";
import { setMobileSearchOpen } from "../../common/store/actions";

const CellGuideMobileHeader = () => {
  const { cellGuideTitle, cellGuideNav, skinnyMode, mobileSearchIsOpen } =
    useContext(StateContext);
  const dispatch = useContext(DispatchContext);

  const [pageNavIsOpen, setPageNavIsOpen] = useState(false);

  const search = (
    <MobileSearchBarWrapper
      onBlur={() => {
        // Hide the input box on blur
        if (dispatch) {
          dispatch(setMobileSearchOpen(false));
        }
      }}
    >
      <CellGuideCardSearchBar autoFocus />
    </MobileSearchBarWrapper>
  );

  const navigation = (
    <>
      {/* Flex Item Left */}
      <div id="cellguide-search-icon">
        <ButtonIcon
          sdsIcon="search"
          id="cellguide-search-icon"
          onClick={() => {
            if (dispatch) {
              dispatch(setMobileSearchOpen(true));
            }
          }}
        />
      </div>

      {/* Flex Item Middle */}
      <StyledTitle id="cellguide-topic">{cellGuideTitle}</StyledTitle>

      {/* Flex Item Right */}
      <div id="cellguide-nav-dropdown">
        <ButtonIcon
          sdsIcon="chevronDown"
          onClick={() => {
            setPageNavIsOpen(!pageNavIsOpen);
          }}
        />
      </div>
    </>
  );

  return (
    <>
      {skinnyMode && (
        <MobileHeaderWrapper>
          {/* CellGuide Header - First Row */}
          <MobileHeader>
            {mobileSearchIsOpen ? search : navigation}
          </MobileHeader>

          {/* CellGuide Page Nav - Second Row */}
          <MobilePageNavWrapper
            onClick={() => {
              // Close nav when clicking on a section
              setPageNavIsOpen(false);
            }}
          >
            {pageNavIsOpen && cellGuideNav}
          </MobilePageNavWrapper>
        </MobileHeaderWrapper>
      )}
    </>
  );
};

export default CellGuideMobileHeader;
