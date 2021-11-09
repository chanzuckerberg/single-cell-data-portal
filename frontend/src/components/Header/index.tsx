import { AnchorButton } from "@blueprintjs/core";
import Link from "next/link";
import { useRouter } from "next/router";
import { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useUserInfo } from "src/common/queries/auth";
import { HomepageLink } from "../common/HomepageLink";
import AuthButtons from "./components/AuthButtons";
import LearnButton from "./components/LearnButton";
import {
  AuthButtonWrapper,
  LearnButtonWrapper,
  Left,
  LinkWrapper,
  MainWrapper,
  Right,
  Wrapper,
} from "./style";

const Header: FC = () => {
  const isCurator = get(FEATURES.CURATOR) === BOOLEAN.TRUE;
  const { data: userInfo } = useUserInfo(isCurator);
  const { asPath } = useRouter();
  const isMyCollectionsActive = asPath === ROUTES.MY_COLLECTIONS;
  const isMyCollectionsShown = userInfo?.name && isCurator;

  return (
    <Wrapper>
      <MainWrapper>
        <Left>
          <HomepageLink />
        </Left>
        <Right>
          {!isMyCollectionsShown && (
            <LinkWrapper>
              <Link href={ROUTES.MY_COLLECTIONS} passHref>
                <AnchorButton
                  active={isMyCollectionsActive}
                  href="passHref"
                  minimal
                  text="My Collections"
                />
              </Link>
            </LinkWrapper>
          )}
          <LearnButtonWrapper>
            <LearnButton />
          </LearnButtonWrapper>
          <AuthButtonWrapper>
            <AuthButtons />
          </AuthButtonWrapper>
        </Right>
      </MainWrapper>
    </Wrapper>
  );
};

export default Header;
