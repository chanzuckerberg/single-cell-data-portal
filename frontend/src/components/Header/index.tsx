import { Intent } from "@blueprintjs/core";
import Link from "next/link";
import { FC } from "react";
import { EXTERNAL_LINKS, ROUTES } from "src/common/constants/routes";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useUserInfo } from "src/common/queries/auth";
import { HomepageLink } from "../common/HomepageLink";
import AuthButtons from "./components/AuthButtons";
import {
  LearnButton,
  MainWrapper,
  MyCollectionsButton,
  Right,
  Wrapper,
} from "./style";

const Header: FC = () => {
  const isCurator = get(FEATURES.CURATOR) === BOOLEAN.TRUE;
  const { data: userInfo } = useUserInfo(isCurator);

  const isMyCollectionsShown = userInfo?.name && isCurator;

  return (
    <Wrapper>
      <MainWrapper>
        <HomepageLink />
        <Right>
          <a
            href={EXTERNAL_LINKS.DOCS_HOME_PAGE}
            target="_blank"
            rel="noopener"
          >
            <LearnButton intent={Intent.PRIMARY} minimal>
              Learn
            </LearnButton>
          </a>
          {isMyCollectionsShown && (
            <Link href={ROUTES.MY_COLLECTIONS} passHref>
              <a href="passHref">
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
