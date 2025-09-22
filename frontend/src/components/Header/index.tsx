import { MenuSelect, Button } from "@czi-sds/components";
import { FC, useState } from "react";
import { HomepageLink } from "../common/HomepageLink";
import AuthButtons from "./components/AuthButtons";
import Nav from "src/components/Header/components/Nav";
import BYODModal from "src/components/BYODModal";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { FEATURES } from "src/common/featureFlags/features";
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

  const [isBYODModalOpen, setIsBYODModalOpen] = useState(false);
  const isBYODEnabled = useFeatureFlag(FEATURES.BYOD);

  const handleBYODClick = () => {
    setIsBYODModalOpen(true);
  };

  const handleBYODClose = () => {
    setIsBYODModalOpen(false);
  };

  return (
    <Wrapper id={HEADER_ID}>
      <MainWrapper>
        <Left>
          <HomepageLink />
          <Nav pathname={pathname} />
        </Left>
        <Right>
          {isBYODEnabled && (
            <Button
              sdsType="primary"
              sdsStyle="square"
              onClick={handleBYODClick}
            >
              Explore Your Data
            </Button>
          )}

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

      <BYODModal open={isBYODModalOpen} onClose={handleBYODClose} />
    </Wrapper>
  );
};

export default Header;
