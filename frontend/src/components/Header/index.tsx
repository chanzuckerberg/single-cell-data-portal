import { Intent } from "@blueprintjs/core";
import Link from "next/link";
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
            <Link href={ROUTES.MY_COLLECTIONS}>
              {/* TODO(thuang): Remove the disable once the issue resolves: https://github.com/vercel/next.js/discussions/8207 */}
              {/* eslint-disable-next-line jsx-a11y/anchor-is-valid */}
              <a>
                <MyCollectionsButton intent={Intent.PRIMARY} minimal>
                  My Collections
                </MyCollectionsButton>
              </a>
            </Link>
          )}
          <AuthButtons />
        </Right>
      </MainWrapper>
    </Wrapper>
  );
};

export default Header;
