import {
  AnchorButton,
  Button,
  Menu,
  MenuDivider,
  MenuItem,
  Popover,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { API } from "src/common/API";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useUserInfo } from "src/common/queries/auth";
import { AuthButtonWrapper } from "src/components/Header/style";
import { API_URL } from "src/configs/configs";
import CuratorAPIKeyGenerator from "../CuratorAPIKeyGenerator";

const AuthButtons = (): JSX.Element | null => {
  const hasAuth = get(FEATURES.CURATOR) === BOOLEAN.TRUE;

  const { data: userInfo, isLoading, error } = useUserInfo(hasAuth);

  if (userInfo && error) {
    // (thuang): Force refresh page to log user out
    window.location.reload();
  }

  if (!hasAuth || isLoading) return null;

  return (
    <AuthButtonWrapper>
      {userInfo?.name ? (
        <LoggedInButtons name={userInfo.name} email={userInfo.email} />
      ) : (
        <LoggedOutButtons />
      )}
    </AuthButtonWrapper>
  );
};

function LoggedInButtons({ name, email }: { name?: string; email?: string }) {
  const authName = isEmail(name) ? "Account" : name;

  return (
    <Popover content={<Content />}>
      <Button minimal rightIcon={IconNames.CHEVRON_DOWN} text={authName} />
    </Popover>
  );

  function Content() {
    const curatorAPIFeature = get(FEATURES.CURATOR_API);
    return (
      <Menu>
        {curatorAPIFeature && (
          <>
            <CuratorAPIKeyGenerator />
            <MenuDivider />
          </>
        )}
        <MenuItem
          data-testid="log-out"
          text="Log Out"
          href={`${API_URL}${API.LOG_OUT}`}
          label={email}
        />
      </Menu>
    );
  }
}

function LoggedOutButtons() {
  return (
    <AnchorButton
      href={`${API_URL}${API.LOG_IN}`}
      minimal
      rightIcon={IconNames.CHEVRON_DOWN}
      text="Log In"
    />
  );
}

function isEmail(name?: string): boolean {
  if (!name) return false;

  return name.includes("@");
}

export default AuthButtons;
