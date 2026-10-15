import { MenuSelect } from "@czi-sds/components";
import { FC } from "react";
import { HomepageLink } from "../common/HomepageLink";
import AuthButtons from "./components/AuthButtons";
import Nav from "src/components/Header/components/Nav";
import { useConnect } from "./connect";
import {
  Left,
  MainWrapper,
  Right,
  StyledInputDropdown,
  StyledPopper,
  Wrapper,
} from "./style";
import { FOOTER_OPTIONS, HEADER_ID } from "./constants";

const Header: FC = () => {
  const {
    pathname,
    dropdownOpen,
    dropdownRef,
    anchorEl,
    handleHelpOpen,
    handleHelpClick,
    handleHelpClose,
  } = useConnect();

  return (
    <Wrapper id={HEADER_ID}>
      <MainWrapper>
        <Left>
          <HomepageLink />
          <Nav pathname={pathname} />
        </Left>
        <Right>
          <StyledInputDropdown
            disabled={false}
            label="Help & Documentation"
            sdsStage="default"
            onClick={handleHelpOpen}
            sdsStyle="minimal"
            sdsType="label"
            data-testid="InputDropdown"
          />

          <StyledPopper
            ref={dropdownRef}
            open={dropdownOpen}
            anchorEl={anchorEl}
          >
            <MenuSelect
              search={false}
              options={FOOTER_OPTIONS}
              onChange={handleHelpClick}
              onClose={handleHelpClose}
            />
          </StyledPopper>

          <AuthButtons />
        </Right>
      </MainWrapper>
    </Wrapper>
  );
};

export default Header;
