import { Intent } from "@blueprintjs/core";
import { Link } from "@reach/router";
import React, { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useUserInfo } from "src/common/queries/auth";
import { HomepageLink } from "../common/HomepageLink";
import AuthButtons from "./components/AuthButtons";
import { MainWrapper, MyCollectionsButton, Right, Wrapper } from "./style";

const Header: FC = () => {
  const isAuth = get(FEATURES.AUTH) === BOOLEAN.TRUE;

  const { data: userInfo } = useUserInfo(isAuth);

  const isMyCollectionsShown =
    userInfo?.name && get(FEATURES.CREATE_COLLECTION) === BOOLEAN.TRUE;

  return (
    <Wrapper>
      <MainWrapper>
        <HomepageLink />
        <Right>
          {isMyCollectionsShown && (
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
