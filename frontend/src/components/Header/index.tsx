import { Intent } from "@blueprintjs/core";
import { Link } from "@reach/router";
import React, { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { HomepageLink } from "../common/HomepageLink";
import AuthButtons from "./components/AuthButtons";
import { MainWrapper, MyCollectionsButton, Right, Wrapper } from "./style";

const Header: FC = () => {
  return (
    <Wrapper>
      <MainWrapper>
        <HomepageLink />
        <Right>
          {get(FEATURES.CREATE_COLLECTION) === BOOLEAN.TRUE && (
            <Link to={ROUTES.MY_COLLECTIONS}>
              <MyCollectionsButton intent={Intent.PRIMARY} minimal>
                My Collections
              </MyCollectionsButton>
            </Link>
          )}
          <AuthButtons />
        </Right>
      </MainWrapper>
    </Wrapper>
  );
};

export default Header;
